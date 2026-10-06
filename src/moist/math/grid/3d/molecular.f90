!> Atom-centered molecular integration grid
!>
!> - `new_molecular_grid`: atomic recipes, partition, pruning, reciprocal settings
!> - `update`: translate cached atomic grids, partition, prune
!> - Reciprocal path: auto-sized NUFFT k-grid, period guard, FINUFFT transform
module moist_math_grid_3d_molecular
   use, intrinsic :: iso_fortran_env, only: output_unit, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_data_atomicrad, only: max_elem
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, &
      & moist_math_grid_atomic_recipe_override_type, element_override_index, &
      & default_molecular_recipe
   use moist_math_grid_atomic_grid, only: moist_math_grid_atomic_type, new_atomic_grid
   use moist_math_grid_3d_partition, only: becke_partition_weights, ssf_partition_weights, &
      & pvoronoi_partition_weights, default_stiffness, default_ssf_a, &
      & default_power_width, partition_becke, partition_ssf, partition_pvoronoi
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type, &
                                      integrand_3d
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
   use finufft_mod, only: finufft_opts
   implicit none(type, external)
   private

   public :: moist_math_grid_3d_molecular_type
   public :: new_molecular_grid
   public :: molecular_grid_set_kgrid
   public :: moist_math_grid_3d_molecular_trafo_type, new_molecular_grid_trafo
   public :: integrand_3d
   public :: default_nufft_tol

   !> Default requested FINUFFT relative tolerance
   real(wp), parameter :: default_nufft_tol = 1.0e-10_wp

   !> Default weight-pruning threshold, bohr^3
   !>
   !> Drop final-grid points with |w| below threshold
   real(wp), parameter :: default_wthr = 1.0e-14_wp

   !> Default real-space spacing of the reciprocal grid, bohr
   real(wp), parameter :: default_dr = 0.5_wp

   !> Two pi, used for the reciprocal-grid spacing dk = 2*pi/(N*dr)
   real(wp), parameter :: two_pi = 8.0_wp*atan(1.0_wp)
   !> (2*pi)^3, the inverse-Fourier-transform normalisation denominator
   real(wp), parameter :: two_pi_cubed = two_pi*two_pi*two_pi
   !> Largest magnitude of a scaled type-1/type-2 coordinate handed to FINUFFT
   !>
   !> Older FINUFFT releases require `[-3pi, 3pi]`
   !> Outside this window, use type 3 for portability
   real(wp), parameter :: type12_coord_limit = 3.0_wp*(0.5_wp*two_pi)
   !> Default extra real-space margin, bohr, on each side of the point box
   !>
   !> Auto-sizing margin to avoid FINUFFT periodic wrap
   real(wp), parameter :: default_kgrid_buffer = 2.0_wp
   !> Safety cap on the auto-sized reciprocal-grid modes per axis
   !>
   !> Far low-weight shells can inflate `N = box/dr` and NUFFT memory
   !> Exceeding cap raises an error
   integer, parameter :: max_kgrid_modes_per_axis = 256

   !> Atom-centered molecular integration grid
   !>
   !> - Owned atomic recipes and molecular partition from `new_molecular_grid`
   !> - Coordinates, count and volume fields inherited from the domain bases
   !> - Per-atom and per-shell results public for reading
   !> - No transform handles; destroy borrowed engines before geometry changes
   type, extends(moist_math_grid_3d_type) :: moist_math_grid_3d_molecular_type
      !> Whether `new_molecular_grid` configured this grid
      logical, private :: constructed = .false.
      !> Atomic recipe of every element without an override
      type(moist_math_grid_atomic_recipe_type), private :: recipe
      !> Per-element recipe overrides, shape (noverride)
      type(moist_math_grid_atomic_recipe_override_type), allocatable, private :: overrides(:)

      !* ---------------------------------- Settings ---------------------------------- *!

      !> Becke polynomial iteration count
      integer, private :: becke_k = default_stiffness
      !> Optional SSF cell parameter, replacing the iterated polynomial
      real(wp), allocatable, private :: ssf_a
      !> Molecular partition scheme
      integer, private :: partition = partition_becke
      !> Power-gap switching half-width in bohr**2
      real(wp), private :: power_width = default_power_width
      !> Quadrature-weight pruning threshold in bohr**3; zero retains all nonzero weights
      real(wp), private :: weight_threshold = default_wthr
      !> Bare partition-weight pruning threshold
      real(wp), private :: pruning_threshold = 0.0_wp
      !> Real-space spacing used to size the reciprocal grid, bohr
      real(wp), private :: dr = default_dr
      !> Period margin on each side of the point cloud, bohr
      real(wp), private :: kbuffer = default_kgrid_buffer
      !> Publish weight-dependent Gaussian widths instead of point potentials
      logical, private :: gaussian = .false.
      !> Width scale, xi0 = xi0_factor / w**(1/3), for positive weights
      real(wp), private :: xi0_factor = 1.0_wp
      !> Configure reciprocal geometry on the first domain update
      logical, private :: auto_kgrid = .true.

      !* --------------------------------- Actual grid -------------------------------- *!

      !> Local atomic grids, one per element met so far, shape (nlocal)
      type(moist_math_grid_atomic_type), allocatable, private :: local(:)
      !> CSR-style atom ownership, shape (nat+1)
      !>
      !> Atom i owns [atom_offset(i), atom_offset(i+1)-1]
      !> First offset 1; final offset ngrid+1
      integer, allocatable :: atom_offset(:)
      !> Per-atom requested radial shell count (before cutoffs), shape (nat)
      integer, allocatable :: nrad_per_atom(:)
      !> Per-atom largest angular point count over its shells, shape (nat)
      integer, allocatable :: nang_per_atom(:)
      !> CSR-style shell ownership, shape (nat+1)
      !>
      !> Atom i shells: [atom_shell_offset(i), atom_shell_offset(i+1)-1]
      !> `shell_*` follows radial order before partition and pruning
      integer, allocatable :: atom_shell_offset(:)
      !> Shell radius (bohr), shape (nshell)
      real(wp), allocatable :: shell_r(:)
      !> Angular point count each shell was built with, shape (nshell)
      !>
      !> Pre-pruning count may exceed retained `atom_offset` range
      integer, allocatable :: shell_npts(:)
      !> Angular exactness degree each shell was built with, shape (nshell)
      integer, allocatable :: shell_degree(:)
      !> Area-based spacing estimate `sqrt(4*pi*r**2/N)` per shell as built (bohr), shape (nshell)
      real(wp), allocatable :: shell_spacing(:)

      !* ------------------ Uniform reciprocal (k) grid for the NUFFT ----------------- *!

      !> Mode counts per axis
      !>
      !> `molecular_grid_set_kgrid` builds uniform enclosing k-box
      !> Before setup: `has_kgrid = false`, `npts_k = 0`; no transforms
      !> `nkx*nky*nkz` modes, per-axis `dk*` spacing
      !> CMCL `modeord = 0`: kx fastest, frequencies -N/2 .. N/2-1
      integer :: nkx = 0, nky = 0, nkz = 0
      !> Reciprocal-space spacings (1/bohr)
      real(wp) :: dkx = 0.0_wp, dky = 0.0_wp, dkz = 0.0_wp
      !> Phase reference: real-space point the transform phases are measured from
      !>
      !> Direct grids use `point(1)`; domain updates use solute centroid
      real(wp) :: kref(3) = 0.0_wp
      !> Requested FINUFFT relative tolerance for transforms on this grid
      !>
      !> Set by constructor or `molecular_grid_set_kgrid`
      !> Transform inherits at `prepare` unless overridden
      !> See `default_nufft_tol` for cost and accuracy
      real(wp), private :: nufft_tol = default_nufft_tol
      !> Symmetric period around an explicit phase reference
      logical, private :: centered_period = .false.
      !> Explicit reference offset from the molecular centroid, bohr
      real(wp), private :: kref_offset(3) = 0.0_wp
      !> Whether the reciprocal grid has been configured
      logical :: has_kgrid = .false.
      !> Latest requested geometry, retained after a period error
      type(structure_type), allocatable, private :: molecule
      !> Realized-geometry generation counter
      !>
      !> Increment on geometry, reciprocal-grid, destroy, or configuration changes
      !> Transform captures generation at `prepare`; `check_ready` rejects changes
      integer, private :: geom_generation = 0
   contains
      procedure :: validate => validate_molecular_grid
      procedure :: update => update_molecular
      procedure :: rebuild => rebuild_molecular
      procedure :: kind_name => molecular_kind_name
      procedure :: has_geometry_dependent_xi0 => molecular_xi0_dependent
      !> Reciprocal-space coordinate of k-point j (1..npts_k; needs k-grid set)
      procedure :: kpoint => molecular_grid_kpoint
      !> Free geometry, per-atom results and cached atomic grids (idempotent)
      procedure :: destroy => molecular_grid_destroy
      !> Print a short summary to the given unit (default output_unit)
      procedure :: info => molecular_grid_info
      !> Allocate the matching transform engine bound to this grid
      procedure :: new_trafo => molecular_grid_new_trafo
      !> Real-space spacing that sizes the reciprocal grid (bohr)
      procedure :: get_dr => molecular_get_dr
      !> Period margin on each side of the point cloud (bohr)
      procedure :: get_kbuffer => molecular_get_kbuffer
      !> Gaussian width scale
      procedure :: get_xi0_factor => molecular_get_xi0_factor
      !> Becke polynomial iteration count
      procedure :: get_becke_k => molecular_get_becke_k
      !> Molecular partition scheme
      procedure :: get_partition => molecular_get_partition
      !> Power-gap half-width
      procedure :: get_power_width => molecular_get_power_width
      !> Quadrature-weight pruning threshold
      procedure :: get_weight_threshold => molecular_get_weight_threshold
      !> Bare partition-weight pruning threshold
      procedure :: get_pruning_threshold => molecular_get_pruning_threshold
      !> Requested FINUFFT relative tolerance
      procedure :: get_nufft_tol => molecular_get_nufft_tol
   end type moist_math_grid_3d_molecular_type

   !> Fourier transform engine for a molecular grid
   !>
   !> FINUFFT NUFFT for nonuniform molecular points
   !> Forward strengths `w_j*f(r_j)` match Cartesian `dV*FFT`
   !> Backward scale `dkx*dky*dkz/(2*pi)^3`
   !> Both use `kref` for Cartesian/`potential_3d` phase convention
   !>
   !> Uniform-box dual lattice permits type 1 forward and type 2 backward
   !> `want_type12 = .false.` selects equivalent, costlier type 3
   !>
   !> `prepare(ntrans)` builds two plans and sorts points once
   !> Plans retain mutable scratch; execute each transform from one thread
   !> Batched calls handle `nv` internally
   !> Worker count fixed at `prepare` from the grid's `team_size`
   !> Reprepare after changing thread count
   type, extends(moist_math_grid_3d_trafo_type) :: moist_math_grid_3d_molecular_trafo_type
      !> Grid this trafo transforms on (not owned; must outlive the trafo)
      class(moist_math_grid_3d_molecular_type), pointer :: grid => null()
      !> Number of simultaneous transforms baked into the plans (0 = unprepared)
      integer :: ntrans = 0
      !> Grid's `geom_generation` at the end of a successful `prepare`
      !>
      !> `check_ready` rejects stale geometry
      integer :: geom_generation = -1
      !> Requested FINUFFT relative tolerance actually used by the plans
      !>
      !> Grid `nufft_tol` unless constructor overrides it
      real(wp) :: nufft_tol = default_nufft_tol
      !> Request the type-1/type-2 fast path (the default)
      !>
      !> `.false.` selects type 3 for elementwise comparison with type 1/2
      logical :: want_type12 = .true.
      !> Whether the prepared plans are actually type 1/2 (`.false.` = type 3)
      !>
      !> Fall back to type 3 outside portable `[-3pi, 3pi]` window
      logical :: is_type12 = .false.
      !> Opaque FINUFFT plan handles, real-space to reciprocal-space (0 = none)
      integer(int64) :: plan_r2k = 0_int64
      !> Opaque FINUFFT plan handle, reciprocal-space to real-space (0 = none)
      integer(int64) :: plan_k2r = 0_int64
      !> Shifted molecular-grid coordinates, length ngrid
      !>
      !> Retained for plan lifetime
      !> Type 3: `xj(j) = xyz(1,j) - kref(1)` in bohr
      !> Type 1/2: the same shifted coordinate scaled by `dkx` into
      !> FINUFFT's radian convention, `xj(j) = (xyz(1,j) - kref(1))*dkx`
      real(c_double), allocatable :: xj(:), yj(:), zj(:)
      !> Reciprocal-grid target/source coordinates, length npts_k
      !>
      !> Type-3 only; type 1/2 uses lattice mode indices
      real(c_double), allocatable :: xk(:), yk(:), zk(:)
      !> Real-space strength/result scratch, shape (ngrid, ntrans), reused every call
      complex(c_double_complex), allocatable :: cj(:, :)
   contains
      procedure :: fft_r2k => molecular_trafo_fft_r2k
      procedure :: fft_k2r => molecular_trafo_fft_k2r
      procedure :: prepare => molecular_trafo_prepare
      procedure :: destroy => molecular_trafo_destroy
   end type moist_math_grid_3d_molecular_trafo_type

contains

   !* ================================================================================= *!
   !*                                   Construction                                    *!
   !* ================================================================================= *!

   !> Configure a molecular grid; geometry is realized only by `update`
   !>
   !> Defaults:
   !>
   !> - `recipe`: `default_molecular_recipe`; midpoint 50, HandyMod to 10 bohr,
   !>   110 Lebedev points per shell
   !> - `becke_k = 3` without SSF, `pruning_threshold = 0` (no partition
   !>   pruning), `weight_threshold = 1e-14` bohr**3, point potentials (`gaussian = .false.`, `xi0_factor = 1`)
   !> - `dr = 0.5` bohr, `kbuffer = 2` bohr, `nufft_tol = 1e-10`,
   !>   `reciprocal = .true.`
   !>
   !> Integration grids: `default_element_recipes` per-element table
   !> Unbounded Becke mapping requires `reciprocal = .false.` or `rcut_upper`
   !>
   !> Copy recipes and overrides; leave grid unconstructed on error
   !> Reconfiguration discards geometry and invalidates bound transform plans
   !>
   !> @param[in,out] self            Configured, geometry-free grid
   !> @param[out] error              Invalid option or incomplete recipe
   !> @param[in]  recipe             Atomic recipe of every element without an override
   !> @param[in]  overrides          Per-element recipe overrides; the first listing an element wins
   !> @param[in]  becke_k            Becke polynomial iteration count (>= 1); not together with `ssf_a`
   !> @param[in]  ssf_a              SSF cell parameter in (0, 1]; replaces the Becke polynomial
   !> @param[in]  pruning_threshold  Bare partition-weight pruning threshold in [0, 1)
   !> @param[in]  dr                 Real-space spacing of the reciprocal grid (bohr, > 0)
   !> @param[in]  kbuffer            Period margin per side (bohr, >= 0)
   !> @param[in]  nufft_tol          Requested FINUFFT relative tolerance (> 0)
   !> @param[in]  gaussian           Publish Gaussian widths `xi0` instead of point potentials
   !> @param[in]  xi0_factor         Width scale, `xi0 = xi0_factor/w**(1/3)` (> 0)
   !> @param[in]  reciprocal         Size the k-grid about the centroid on the first update
   !> @param[in]  partition          Scheme constant; default Becke, or SSF when ssf_a is supplied
   !> @param[in]  power_width        Positive power-gap half-width in bohr**2; only with partition_pvoronoi
   !> @param[in]  weight_threshold   Nonnegative quadrature-weight cutoff in bohr**3; zero for exact-zero pruning
   subroutine new_molecular_grid(self, error, recipe, overrides, becke_k, ssf_a, pruning_threshold, &
         & dr, kbuffer, nufft_tol, gaussian, xi0_factor, reciprocal, partition, power_width, weight_threshold)
      !> Configured grid
      type(moist_math_grid_3d_molecular_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Atomic recipe of every element without an override
      type(moist_math_grid_atomic_recipe_type), intent(in), optional :: recipe
      !> Per-element recipe overrides
      type(moist_math_grid_atomic_recipe_override_type), intent(in), optional :: overrides(:)
      !> Becke polynomial iteration count
      integer, intent(in), optional :: becke_k
      !> SSF cell parameter
      real(wp), intent(in), optional :: ssf_a
      !> Bare partition-weight pruning threshold
      real(wp), intent(in), optional :: pruning_threshold
      !> Real-space spacing of the reciprocal grid
      real(wp), intent(in), optional :: dr
      !> Period margin per side
      real(wp), intent(in), optional :: kbuffer
      !> Requested FINUFFT relative tolerance
      real(wp), intent(in), optional :: nufft_tol
      !> Publish Gaussian widths
      logical, intent(in), optional :: gaussian
      !> Gaussian width scale
      real(wp), intent(in), optional :: xi0_factor
      !> Size the reciprocal grid on the first update
      logical, intent(in), optional :: reciprocal

      !> Molecular partition scheme
      integer, intent(in), optional :: partition
      !> Power-gap half-width
      real(wp), intent(in), optional :: power_width
      !> Quadrature-weight pruning threshold
      real(wp), intent(in), optional :: weight_threshold

      integer :: generation

      ! Preserve counter across reset to invalidate earlier transform plans
      generation = self%geom_generation
      call reset_molecular_grid(self)
      self%geom_generation = generation + 1

      if (present(becke_k) .and. present(ssf_a)) then
         call fatal_error(error, "molecular domain: give either becke_k or ssf_a, not both")
         return
      end if

      if (present(ssf_a)) self%partition = partition_ssf
      if (present(partition)) self%partition = partition
      if (present(becke_k) .and. self%partition /= partition_becke) then
         call fatal_error(error, "molecular domain: becke_k requires the Becke partition")
         return
      end if
      if (present(ssf_a) .and. self%partition /= partition_ssf) then
         call fatal_error(error, "molecular domain: ssf_a requires the SSF partition")
         return
      end if
      if (present(power_width) .and. self%partition /= partition_pvoronoi) then
         call fatal_error(error, "molecular domain: power_width requires the power-Voronoi partition")
         return
      end if

      if (present(recipe)) then
         self%recipe = recipe
      else
         call default_molecular_recipe(self%recipe, error)
         if (allocated(error)) return
      end if
      if (present(overrides)) self%overrides = overrides
      if (present(becke_k)) self%becke_k = becke_k
      if (self%partition == partition_ssf) self%ssf_a = default_ssf_a
      if (present(ssf_a)) self%ssf_a = ssf_a
      if (present(power_width)) self%power_width = power_width
      if (present(weight_threshold)) self%weight_threshold = weight_threshold
      if (present(pruning_threshold)) self%pruning_threshold = pruning_threshold
      if (present(dr)) self%dr = dr
      if (present(kbuffer)) self%kbuffer = kbuffer
      if (present(nufft_tol)) self%nufft_tol = nufft_tol
      if (present(gaussian)) self%gaussian = gaussian
      if (present(xi0_factor)) self%xi0_factor = xi0_factor
      if (present(reciprocal)) self%auto_kgrid = reciprocal

      call check_molecular_settings(self, error)
      if (allocated(error)) return
      self%constructed = .true.
   end subroutine new_molecular_grid

   !> Reset a grid to defaults
   !>
   !> @param[out] self Grid, reset by its `intent(out)`
   subroutine reset_molecular_grid(self)
      !> Grid to reset
      type(moist_math_grid_3d_molecular_type), intent(out) :: self
   end subroutine reset_molecular_grid

   !> Validate the stored configuration
   !>
   !> @param[in] self Grid configuration
   !> @param[out] error Unconstructed grid or invalid setting
   subroutine validate_molecular_grid(self, error)
      !> Grid configuration
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Unconstructed grid or invalid setting
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%constructed) then
         call fatal_error(error, "molecular domain: construct the grid with new_molecular_grid first")
         return
      end if
      call check_molecular_settings(self, error)
   end subroutine validate_molecular_grid

   !> Validate options and recipes before grid construction
   !>
   !> @param[in] self Grid configuration
   !> @param[out] error Invalid setting
   subroutine check_molecular_settings(self, error)
      !> Grid configuration
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Invalid setting
      type(error_type), allocatable, intent(out) :: error

      character(len=32) :: label
      integer :: i

      if (.not. all(ieee_is_finite([self%pruning_threshold, self%dr, self%kbuffer, self%nufft_tol, &
                                    self%xi0_factor, self%power_width, self%weight_threshold]))) then
         call fatal_error(error, "molecular domain: settings must be finite")
      else if (.not. any(self%partition == [partition_becke, partition_ssf, partition_pvoronoi])) then
         call fatal_error(error, "molecular domain: unknown partition scheme")
      else if (self%power_width <= 0.0_wp) then
         call fatal_error(error, "molecular domain: power_width must be positive")
      else if (self%weight_threshold < 0.0_wp) then
         call fatal_error(error, "molecular domain: weight_threshold must be nonnegative")
      else if (self%becke_k < 1) then
         call fatal_error(error, "molecular domain: becke_k must be >= 1")
      else if (self%pruning_threshold < 0.0_wp .or. self%pruning_threshold >= 1.0_wp) then
         call fatal_error(error, "molecular domain: pruning_threshold must be in [0, 1)")
      else if (self%dr <= 0.0_wp .or. self%kbuffer < 0.0_wp .or. self%nufft_tol <= 0.0_wp) then
         call fatal_error(error, "molecular domain: invalid dr, kbuffer or nufft_tol")
      else if (self%xi0_factor <= 0.0_wp) then
         call fatal_error(error, "molecular domain: xi0_factor must be positive")
      end if
      if (allocated(error)) return
      if (allocated(self%ssf_a)) then
         if (.not. ieee_is_finite(self%ssf_a)) then
            call fatal_error(error, "molecular domain: ssf_a must be finite")
         else if (self%ssf_a <= 0.0_wp .or. self%ssf_a > 1.0_wp) then
            call fatal_error(error, "molecular domain: ssf_a must be in (0, 1]")
         end if
      end if
      if (allocated(error)) return

      call check_atomic_recipe(self%recipe, "default recipe", error)
      if (allocated(error)) return
      if (.not. allocated(self%overrides)) return
      do i = 1, size(self%overrides)
         write (label, "(a,i0)") "override ", i
         if (.not. allocated(self%overrides(i)%elements)) then
            call fatal_error(error, "molecular domain: "//trim(label)//" lists no elements")
            return
         end if
         if (size(self%overrides(i)%elements) < 1) then
            call fatal_error(error, "molecular domain: "//trim(label)//" lists no elements")
            return
         end if
         call check_atomic_recipe(self%overrides(i)%recipe, trim(label), error)
         if (allocated(error)) return
      end do
   end subroutine check_molecular_settings

   !> Check atomic recipe completeness and radial settings
   !>
   !> Check rule, mapping, and generator when first building each element
   !>
   !> @param[in]  recipe  Atomic recipe
   !> @param[in]  label   Recipe name for the message
   !> @param[out] error   Incomplete recipe, bad node count, cutoffs, or shell policy
   subroutine check_atomic_recipe(recipe, label, error)
      !> Atomic recipe
      type(moist_math_grid_atomic_recipe_type), intent(in) :: recipe
      !> Recipe name for the message
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(error_type), allocatable :: shell_error

      if (.not. (allocated(recipe%radial%rule) .and. allocated(recipe%radial%mapping) &
         & .and. allocated(recipe%angular) .and. allocated(recipe%shells))) then
         call fatal_error(error, "molecular domain: "//label//" needs a radial rule, a radial "// &
            & "mapping, an angular generator and a shell policy")
         return
      end if
      if (recipe%radial%npts < 1) then
         call fatal_error(error, "molecular domain: "//label//" needs at least one radial node")
         return
      end if
      if (allocated(recipe%radial%rcut_lower)) then
         if (.not. ieee_is_finite(recipe%radial%rcut_lower)) then
            call fatal_error(error, "molecular domain: settings must be finite ("//label//" rcut_lower)")
            return
         end if
      end if
      if (allocated(recipe%radial%rcut_upper)) then
         if (.not. ieee_is_finite(recipe%radial%rcut_upper)) then
            call fatal_error(error, "molecular domain: settings must be finite ("//label//" rcut_upper)")
            return
         end if
      end if
      if (allocated(recipe%radial%rcut_lower) .and. allocated(recipe%radial%rcut_upper)) then
         if (.not. recipe%radial%rcut_lower < recipe%radial%rcut_upper) then
            call fatal_error(error, "molecular domain: "//label//" needs rcut_lower < rcut_upper")
            return
         end if
      end if
      call recipe%shells%validate(shell_error)
      if (allocated(shell_error)) then
         call fatal_error(error, "molecular domain: "//label//": "//shell_error%message)
         return
      end if
   end subroutine check_atomic_recipe

   !* ================================================================================= *!
   !*                                 Geometry updates                                  *!
   !* ================================================================================= *!

   !> Commit a constructed molecular discretization
   !>
   !> @param[in,out] source Candidate grid
   !> @param[in,out] dest Destination retaining settings, molecule and context
   subroutine move_molecular_geometry(source, dest)
      !> Candidate grid
      type(moist_math_grid_3d_molecular_type), intent(inout) :: source
      !> Destination grid
      class(moist_math_grid_3d_molecular_type), intent(inout) :: dest
      dest%ngrid = source%ngrid
      dest%natom = source%natom
      dest%npts_k = source%npts_k
      dest%nkx = source%nkx
      dest%nky = source%nky
      dest%nkz = source%nkz
      dest%dkx = source%dkx
      dest%dky = source%dky
      dest%dkz = source%dkz
      dest%kref = source%kref
      dest%has_kgrid = source%has_kgrid
      dest%centered_period = source%centered_period
      dest%kref_offset = source%kref_offset
      ! New point cloud invalidates earlier transform plans
      dest%geom_generation = dest%geom_generation + 1
      call move_alloc(source%xyz, dest%xyz)
      call move_alloc(source%w, dest%w)
      call move_alloc(source%owner, dest%owner)
      call move_alloc(source%xi0, dest%xi0)
      call move_alloc(source%atom_offset, dest%atom_offset)
      call move_alloc(source%nrad_per_atom, dest%nrad_per_atom)
      call move_alloc(source%nang_per_atom, dest%nang_per_atom)
      call move_alloc(source%atom_shell_offset, dest%atom_shell_offset)
      call move_alloc(source%shell_r, dest%shell_r)
      call move_alloc(source%shell_npts, dest%shell_npts)
      call move_alloc(source%shell_degree, dest%shell_degree)
      call move_alloc(source%shell_spacing, dest%shell_spacing)
   end subroutine move_molecular_geometry

   !> Move atomic grids, recompute partition weights and guard the fixed period
   !>
   !> @param[in,out] self Domain instance
   !> @param[in] mol Solute structure
   !> @param[out] error Unconstructed grid, invalid geometry or insufficient period
   subroutine update_molecular(self, mol, error)
      !> Domain instance
      class(moist_math_grid_3d_molecular_type), intent(inout) :: self
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Unconstructed grid, invalid geometry or insufficient period
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%constructed) then
         call fatal_error(error, "molecular domain: construct the grid with new_molecular_grid before update")
         return
      end if
      if (mol%nat < 1) then
         call fatal_error(error, "molecular domain: at least one solute atom is required")
         return
      end if
      if (.not. all(ieee_is_finite(mol%xyz))) then
         call fatal_error(error, "molecular domain: solute coordinates must be finite")
         return
      end if
      ! Partition requires each atom's covalent radius
      if (any(mol%num < 1) .or. any(mol%num > max_elem)) then
         call fatal_error(error, "molecular domain: atomic numbers must be between 1 and the "// &
            & "last tabulated element")
         return
      end if
      self%molecule = mol
      call realize_molecular(self, .false., error)
   end subroutine update_molecular

   !> Rebuild the reciprocal period for the latest requested geometry
   !>
   !> @param[in,out] self Domain instance
   !> @param[out] error Invalid grid or no previous update
   subroutine rebuild_molecular(self, error)
      !> Domain instance
      class(moist_math_grid_3d_molecular_type), intent(inout) :: self
      !> Invalid grid or no previous update
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(self%molecule)) then
         call fatal_error(error, "molecular domain: update before rebuild")
         return
      end if
      call realize_molecular(self, .true., error)
   end subroutine rebuild_molecular

   !> Domain name for diagnostics, including the inherited gradient stub
   !>
   !> @param[in] self Domain instance
   function molecular_kind_name(self) result(name)
      !> Domain instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Diagnostic name
      character(len=:), allocatable :: name
      name = "molecular"
   end function molecular_kind_name

   !> Construct candidate geometry, committing only after all guards pass
   !>
   !> Cached local grids preserve quadrature points across geometries
   !> Reapplied threshold pruning may change point counts
   !>
   !> @param[in,out] self Domain with a requested molecular geometry
   !> @param[in] reset_period Allow the reciprocal grid to change
   !> @param[out] error Construction or period error
   subroutine realize_molecular(self, reset_period, error)
      !> Domain with a requested geometry
      class(moist_math_grid_3d_molecular_type), intent(inout) :: self
      !> Allow the reciprocal grid to change
      logical, intent(in) :: reset_period
      !> Construction or period error
      type(error_type), allocatable, intent(out) :: error
      !> Candidate point cloud
      type(moist_math_grid_3d_molecular_type), allocatable :: grid
      !> Centroid and required symmetric period
      real(wp) :: center(3), reference(3), required(3), period(3), tolerance(3)
      !> Allocation status
      integer :: stat
      !> Whether to retain the previous reciprocal grid
      logical :: preserve_period

      call self%validate(error)
      if (allocated(error)) return
      allocate (grid, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "molecular domain: cannot allocate grid")
         return
      end if
      center = sum(self%molecule%xyz, dim=2)/real(self%molecule%nat, wp)
      call assemble_molecular_grid(self, grid, error)
      if (allocated(error)) return
      if (grid%ngrid < 1) then
         call fatal_error(error, "molecular domain: pruning removed every point")
         return
      end if
      if (.not. all(ieee_is_finite(grid%w)) .or. .not. all(ieee_is_finite(grid%xyz))) then
         call fatal_error(error, "molecular domain: quadrature produced non-finite points or weights")
         return
      end if
      if (self%gaussian) then
         if (any(grid%w <= 0.0_wp)) then
            call fatal_error(error, "molecular domain: Gaussian widths require positive integration weights")
            return
         end if
      end if
      preserve_period = self%has_kgrid .and. .not. reset_period
      if (preserve_period) then
         if (self%centered_period) then
            reference = center + self%kref_offset
            required = 2.0_wp*(max(abs(maxval(grid%xyz, dim=2) - reference), &
               & abs(minval(grid%xyz, dim=2) - reference)) + self%kbuffer)
         else
            reference = grid%xyz(:, 1)
            required = maxval(grid%xyz, dim=2) - minval(grid%xyz, dim=2) + 2.0_wp*self%kbuffer
         end if
         period = two_pi/[self%dkx, self%dky, self%dkz]
         tolerance = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp, period)
         if (any(required > period + tolerance)) then
            call fatal_error(error, "molecular domain: fixed reciprocal period exceeded; call rebuild")
            return
         end if
         grid%nkx = self%nkx
         grid%nky = self%nky
         grid%nkz = self%nkz
         grid%dkx = self%dkx
         grid%dky = self%dky
         grid%dkz = self%dkz
         grid%npts_k = self%npts_k
         grid%nufft_tol = self%nufft_tol
         grid%kref = reference
         grid%centered_period = self%centered_period
         grid%kref_offset = self%kref_offset
         grid%has_kgrid = .true.
      else if (self%auto_kgrid .or. (reset_period .and. self%has_kgrid)) then
         if (self%centered_period .or. self%auto_kgrid) then
            call molecular_grid_set_kgrid(grid, self%dr, error, buffer=self%kbuffer, &
               & nufft_tol=self%nufft_tol, reference=center + self%kref_offset)
         else
            call molecular_grid_set_kgrid(grid, self%dr, error, buffer=self%kbuffer, &
               & nufft_tol=self%nufft_tol)
         end if
         if (allocated(error)) return
      end if

      if (self%gaussian) then
         allocate (grid%xi0(grid%ngrid), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "molecular domain: cannot allocate Gaussian widths")
            return
         end if
         grid%xi0 = self%xi0_factor/grid%w**(1.0_wp/3.0_wp)
      end if
      call move_molecular_geometry(grid, self)
   end subroutine realize_molecular

   !> Assemble a candidate point cloud
   !>
   !> - Atom order, then radial-shell order, then angular-table order
   !> - Point coordinates `r*u + R_A` per component
   !> - Weight: local `r^2 dr dOmega` times owner partition weight
   !> - Prune `|w| < 1e-14` and bare weight below positive threshold
   !>
   !> @param[in,out] self   Configured grid with a requested geometry; caches local grids
   !> @param[in,out] grid   Candidate receiving points, weights, and per-atom results
   !> @param[out]    error  Atomic grid, point count, or allocation failure
   subroutine assemble_molecular_grid(self, grid, error)
      !> Configured grid with a requested geometry
      class(moist_math_grid_3d_molecular_type), intent(inout) :: self
      !> Candidate grid
      type(moist_math_grid_3d_molecular_type), intent(inout) :: grid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Atomic number per atom
      integer, allocatable :: numbers(:)
      !> Cache slot of each atom's local grid
      integer, allocatable :: slot(:)
      !> Requested shell count and largest angular count per atom
      integer, allocatable :: nrad_atom(:), nang_atom(:)
      !> Raw (unpruned) points, weights, partition weights, and owners
      real(wp), allocatable :: raw_xyz(:, :), raw_w(:), raw_bw(:)
      integer, allocatable :: raw_atom(:)
      integer(int64) :: total
      integer :: nat, iat, i, ig, lo, nraw, nshell, ish, stat
      real(wp) :: r

      nat = self%molecule%nat
      allocate (numbers(nat), slot(nat), nrad_atom(nat), nang_atom(nat))
      do iat = 1, nat
         numbers(iat) = self%molecule%num(self%molecule%id(iat))
         call local_grid_slot(self, numbers(iat), slot(iat), error)
         if (allocated(error)) return
      end do

      total = 0_int64
      nshell = 0
      do iat = 1, nat
         associate (local => self%local(slot(iat)))
            total = total + int(local%npts, int64)
            nshell = nshell + local%nshell
            nrad_atom(iat) = local%nshell_requested
            nang_atom(iat) = maxval(local%shell_npts)
         end associate
      end do
      call check_point_count(total, error)
      if (allocated(error)) return
      nraw = int(total)

      allocate (raw_xyz(3, nraw), stat=stat)
      if (stat == 0) allocate (raw_w(nraw), stat=stat)
      if (stat == 0) allocate (raw_bw(nraw), stat=stat)
      if (stat == 0) allocate (raw_atom(nraw), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "molecular domain: cannot allocate the raw point buffer")
         return
      end if

      ig = 0
      do iat = 1, nat
         associate (local => self%local(slot(iat)), mol => self%molecule)
            lo = ig + 1
            do i = 1, local%npts
               ig = ig + 1
               r = local%shell_r(local%shell(i))
               raw_xyz(1, ig) = r*local%u(1, i) + mol%xyz(1, iat)
               raw_xyz(2, ig) = r*local%u(2, i) + mol%xyz(2, iat)
               raw_xyz(3, ig) = r*local%u(3, i) + mol%xyz(3, iat)
               raw_atom(ig) = iat
            end do
            if (ig >= lo) then
               ! Parallel points; serial atom loop
               select case (self%partition)
               case (partition_becke)
                  call becke_partition_weights(iat, raw_xyz(:, lo:ig), mol%xyz, numbers, raw_bw(lo:ig), &
                     & stiffness=self%becke_k, nthreads=self%team_size())
               case (partition_ssf)
                  call ssf_partition_weights(iat, raw_xyz(:, lo:ig), mol%xyz, numbers, raw_bw(lo:ig), a=self%ssf_a, &
                     & nthreads=self%team_size())
               case (partition_pvoronoi)
                  call pvoronoi_partition_weights(iat, raw_xyz(:, lo:ig), mol%xyz, numbers, raw_bw(lo:ig), &
                     & width=self%power_width, nthreads=self%team_size())
               case default
                  call fatal_error(error, "molecular domain: unknown partition scheme")
                  return
               end select
            end if
            do i = 1, local%npts
               raw_w(lo + i - 1) = local%w(i)*raw_bw(lo + i - 1)
            end do
         end associate
      end do

      call finalise_grid(grid, nat, nrad_atom, nang_atom, nraw, raw_xyz, raw_w, raw_atom, &
         & self%weight_threshold, raw_bw=raw_bw, bthr=self%pruning_threshold)
      deallocate (raw_xyz, raw_w, raw_bw, raw_atom)

      allocate (grid%atom_shell_offset(nat + 1), stat=stat)
      if (stat == 0) allocate (grid%shell_r(nshell), stat=stat)
      if (stat == 0) allocate (grid%shell_npts(nshell), stat=stat)
      if (stat == 0) allocate (grid%shell_degree(nshell), stat=stat)
      if (stat == 0) allocate (grid%shell_spacing(nshell), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "molecular domain: cannot allocate the per-shell results")
         return
      end if
      ish = 0
      do iat = 1, nat
         associate (local => self%local(slot(iat)))
            grid%atom_shell_offset(iat) = ish + 1
            grid%shell_r(ish + 1:ish + local%nshell) = local%shell_r
            grid%shell_npts(ish + 1:ish + local%nshell) = local%shell_npts
            grid%shell_degree(ish + 1:ish + local%nshell) = local%shell_degree
            grid%shell_spacing(ish + 1:ish + local%nshell) = local%shell_spacing
            ish = ish + local%nshell
         end associate
      end do
      grid%atom_shell_offset(nat + 1) = ish + 1
      ! Measure candidate reciprocal offset from current centroid
      grid%molecule = self%molecule
   end subroutine assemble_molecular_grid

   !> Cache slot of the local atomic grid of element z, built on first use
   !>
   !> First override for z, otherwise default recipe
   !>
   !> @param[in,out] self   Configured grid owning the cache
   !> @param[in]     z      Atomic number
   !> @param[out]    slot   Index into the cache
   !> @param[out]    error  Atomic grid construction failure
   subroutine local_grid_slot(self, z, slot, error)
      !> Configured grid owning the cache
      class(moist_math_grid_3d_molecular_type), intent(inout) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Index into the cache
      integer, intent(out) :: slot
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_type) :: atom
      type(moist_math_grid_atomic_type), allocatable :: grown(:)
      type(error_type), allocatable :: atom_error
      integer :: k, n, idx

      slot = 0
      n = 0
      if (allocated(self%local)) then
         n = size(self%local)
         do k = 1, n
            if (self%local(k)%z == z) then
               slot = k
               return
            end if
         end do
      end if

      idx = element_override_index(z, self%overrides)
      if (idx == 0) then
         call new_atomic_grid(atom, self%recipe, z, atom_error)
      else
         call new_atomic_grid(atom, self%overrides(idx)%recipe, z, atom_error)
      end if
      if (allocated(atom_error)) then
         call fatal_error(error, "molecular domain: "//atom_error%message)
         return
      end if

      allocate (grown(n + 1))
      do k = 1, n
         grown(k) = self%local(k)
      end do
      grown(n + 1) = atom
      call move_alloc(grown, self%local)
      slot = n + 1
   end subroutine local_grid_slot

   !> Report the weight dependence of optional Gaussian widths
   !>
   !> @param[in] self Domain instance
   function molecular_xi0_dependent(self) result(dependent)
      !> Domain instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Whether gradient-phase potentials need dphi_dxi
      logical :: dependent
      dependent = self%gaussian
   end function molecular_xi0_dependent

   !> Guard a raw-buffer point-count accumulator against integer overflow
   !>
   !> Caller sums atom counts in int64 before default-integer allocation
   !>
   !> @param[in]  total  Accumulated point count
   !> @param[out] error  Set when `total` exceeds `huge(0)`
   subroutine check_point_count(total, error)
      !> Accumulated point count
      integer(int64), intent(in) :: total
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (total > int(huge(0), int64)) then
         call fatal_error(error, "molecular grid: point count exceeds the representable integer range")
      end if
   end subroutine check_point_count

   !> Compact the raw grid buffer into the candidate grid
   !>
   !> Drop zero and `|w| < wthr` points; rebuild `atom_offset`
   !>
   !> @param[out] self       Grid to fill
   !> @param[in]  nat        Number of atoms
   !> @param[in]  nrad_atom  Per-atom requested radial sizes, shape (nat)
   !> @param[in]  nang_atom  Per-atom largest angular point counts, shape (nat)
   !> @param[in]  nraw       Number of raw (pre-pruning) points
   !> @param[in]  raw_xyz    Raw point coordinates, shape (3, nraw), bohr
   !> @param[in]  raw_w      Raw quadrature weights, shape (nraw), bohr^3
   !> @param[in]  raw_atom   Raw point owning atom, shape (nraw), one based
   !> @param[in]  wthr       Pruning threshold on |raw_w| (bohr^3)
   !> @param[in]  raw_bw     Optional bare partition weight of each raw point, shape (nraw)
   !> @param[in]  bthr       Optional partition-weight prune threshold
   subroutine finalise_grid(self, nat, nrad_atom, nang_atom, &
         & nraw, raw_xyz, raw_w, raw_atom, wthr, raw_bw, bthr)
      !> Grid to fill
      type(moist_math_grid_3d_molecular_type), intent(out) :: self
      !> Number of atoms
      integer, intent(in) :: nat
      !> Per-atom requested radial sizes, shape (nat)
      integer, intent(in) :: nrad_atom(:)
      !> Per-atom largest angular point counts, shape (nat)
      integer, intent(in) :: nang_atom(:)
      !> Number of raw (pre-pruning) points
      integer, intent(in) :: nraw
      !> Raw point coordinates, shape (3, nraw), bohr
      real(wp), intent(in) :: raw_xyz(:, :)
      !> Raw quadrature weights, shape (nraw), bohr^3
      real(wp), intent(in) :: raw_w(:)
      !> Raw point owning atom, shape (nraw), one based
      integer, intent(in) :: raw_atom(:)
      !> Pruning threshold on |raw_w| (bohr^3)
      real(wp), intent(in) :: wthr
      !> Bare partition weight of each raw point, shape (nraw)
      !>
      !> Dimensionless `w_A(r)` before radial Jacobian
      !> Sums to 1 across atoms; suitable near nuclei
      real(wp), intent(in), optional :: raw_bw(:)
      !> Partition-weight prune threshold
      !>
      !> Drop points owned more strongly by another atom
      !> Zero or absent adds no pruning beyond `wthr`
      real(wp), intent(in), optional :: bthr

      integer :: i, iat, ngrid
      logical :: prune_bw
      logical, allocatable :: keep(:)

      ! Explicit threshold opt-in; absent threshold never prunes by bare weight
      prune_bw = .false.
      if (present(raw_bw) .and. present(bthr)) then
         if (bthr > 0.0_wp) prune_bw = .true.
      end if

      ! Retain NaN weights for caller's finiteness check
      allocate (keep(nraw))
      if (prune_bw) then
         do i = 1, nraw
            keep(i) = .not. (raw_w(i) == 0.0_wp .or. abs(raw_w(i)) < wthr) .and. .not. (raw_bw(i) < bthr)
         end do
      else
         do i = 1, nraw
            keep(i) = .not. (raw_w(i) == 0.0_wp .or. abs(raw_w(i)) < wthr)
         end do
      end if

      ! First pass: count retained points
      ngrid = 0
      do i = 1, nraw
         if (keep(i)) ngrid = ngrid + 1
      end do

      self%ngrid = ngrid
      self%natom = nat
      allocate (self%xyz(3, ngrid), self%w(ngrid), self%owner(ngrid))
      allocate (self%atom_offset(nat + 1))
      allocate (self%nrad_per_atom(nat), self%nang_per_atom(nat))

      self%nrad_per_atom = nrad_atom
      self%nang_per_atom = nang_atom

      ! Second pass: copy retained points and build atom_offset
      ! `raw_atom` increases by atom order
      ! End sentinel gives empty ranges to atoms without points
      self%atom_offset = ngrid + 1
      self%atom_offset(1) = 1
      ngrid = 0
      iat = 1
      do i = 1, nraw
         ! Advance atom_offset for any atoms that start at this position
         do while (iat < raw_atom(i))
            iat = iat + 1
            self%atom_offset(iat) = ngrid + 1
         end do
         if (keep(i)) then
            ngrid = ngrid + 1
            self%xyz(:, ngrid) = raw_xyz(:, i)
            self%w(ngrid) = raw_w(i)
            self%owner(ngrid) = raw_atom(i)
         end if
      end do
      ! Fill remaining (possibly empty) tail atoms
      do while (iat < nat)
         iat = iat + 1
         self%atom_offset(iat) = ngrid + 1
      end do
      self%atom_offset(nat + 1) = ngrid + 1
   end subroutine finalise_grid

   !* ================================================================================= *!
   !*                                Settings accessors                                 *!
   !* ================================================================================= *!

   !> Real-space spacing that sizes the reciprocal grid
   !>
   !> @param[in] self Grid instance
   pure function molecular_get_dr(self) result(dr)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Spacing (bohr)
      real(wp) :: dr
      dr = self%dr
   end function molecular_get_dr

   !> Period margin on each side of the point cloud
   !>
   !> @param[in] self Grid instance
   pure function molecular_get_kbuffer(self) result(kbuffer)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Margin (bohr)
      real(wp) :: kbuffer
      kbuffer = self%kbuffer
   end function molecular_get_kbuffer

   !> Gaussian width scale, `xi0 = xi0_factor/w**(1/3)`
   !>
   !> @param[in] self Grid instance
   pure function molecular_get_xi0_factor(self) result(xi0_factor)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Width scale
      real(wp) :: xi0_factor
      xi0_factor = self%xi0_factor
   end function molecular_get_xi0_factor

   !> Becke polynomial iteration count
   !>
   !> @param[in] self Grid instance
   pure function molecular_get_becke_k(self) result(becke_k)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Iteration count
      integer :: becke_k
      becke_k = self%becke_k
   end function molecular_get_becke_k

   !> Molecular partition scheme
   !>
   !> @param[in] self Grid instance
   pure function molecular_get_partition(self) result(value)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Molecular partition scheme
      integer :: value
      value = self%partition
   end function molecular_get_partition

   !> Power-gap switching half-width in bohr**2
   !>
   !> @param[in] self Grid instance
   pure function molecular_get_power_width(self) result(value)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Power-gap switching half-width in bohr**2
      real(wp) :: value
      value = self%power_width
   end function molecular_get_power_width

   !> Quadrature-weight pruning threshold in bohr**3
   !>
   !> @param[in] self Grid instance
   pure function molecular_get_weight_threshold(self) result(value)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Quadrature-weight pruning threshold in bohr**3
      real(wp) :: value
      value = self%weight_threshold
   end function molecular_get_weight_threshold

   !> Bare partition-weight pruning threshold
   !>
   !> @param[in] self Grid instance
   pure function molecular_get_pruning_threshold(self) result(pruning_threshold)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Threshold in [0, 1)
      real(wp) :: pruning_threshold
      pruning_threshold = self%pruning_threshold
   end function molecular_get_pruning_threshold

   !> Requested FINUFFT relative tolerance
   !>
   !> @param[in] self Grid instance
   pure function molecular_get_nufft_tol(self) result(nufft_tol)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Relative tolerance
      real(wp) :: nufft_tol
      nufft_tol = self%nufft_tol
   end function molecular_get_nufft_tol

   !* ================================================================================= *!
   !*                                      Destroy                                      *!
   !* ================================================================================= *!

   !> Free geometry, per-atom results, and atomic caches
   !>
   !> Keep configuration for later `update`
   !>
   !> @param[in,out] self  Grid instance
   subroutine molecular_grid_destroy(self)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(inout) :: self

      call self%clear_geometry()
      if (allocated(self%molecule)) deallocate (self%molecule)
      if (allocated(self%atom_offset)) deallocate (self%atom_offset)
      if (allocated(self%nrad_per_atom)) deallocate (self%nrad_per_atom)
      if (allocated(self%nang_per_atom)) deallocate (self%nang_per_atom)
      if (allocated(self%atom_shell_offset)) deallocate (self%atom_shell_offset)
      if (allocated(self%shell_r)) deallocate (self%shell_r)
      if (allocated(self%shell_npts)) deallocate (self%shell_npts)
      if (allocated(self%shell_degree)) deallocate (self%shell_degree)
      if (allocated(self%shell_spacing)) deallocate (self%shell_spacing)
      if (allocated(self%local)) deallocate (self%local)
      self%ngrid = 0
      self%npts_k = 0
      self%nkx = 0; self%nky = 0; self%nkz = 0
      self%dkx = 0.0_wp; self%dky = 0.0_wp; self%dkz = 0.0_wp
      self%kref = 0.0_wp
      self%has_kgrid = .false.
      self%centered_period = .false.
      self%kref_offset = 0.0_wp
      ! Invalidate earlier bound transform plans
      self%geom_generation = self%geom_generation + 1
   end subroutine molecular_grid_destroy

   !* ================================================================================= *!
   !*                             Reciprocal grid and trafo                             *!
   !* ================================================================================= *!

   !> Configure the uniform reciprocal (k) grid used by the NUFFT
   !>
   !> - `dr` fixes Nyquist reach `k_max = pi/dr`; box encloses all points
   !> - Per axis: `L = extent + 2*buffer`; smallest even N with `N*dr >= L`
   !> - Reciprocal spacing `dk = 2*pi/(N*dr)`
   !> - CMCL modes: kx fastest, frequencies -N/2 .. N/2-1
   !> - Default phase reference `point(1)` matches Cartesian/`potential_3d`
   !> - Automatic reciprocal setup uses solute centroid
   !> - Explicit setup preserves chosen reference across updates
   !>
   !> @param[in,out] self       Grid whose real-space points are already built
   !> @param[in]     dr         Real-space spacing of the implied box (bohr, > 0)
   !> @param[out]    error      Propagated error (no points, dr <= 0, ...)
   !> @param[in]     buffer     Optional margin per side (bohr); default 2 bohr
   !> @param[in]     nufft_tol  Optional requested FINUFFT relative tolerance;
   !>                           default `default_nufft_tol` (1e-10)
   !> @param[in]     reference  Optional phase origin; sizes a symmetric period about it
   subroutine molecular_grid_set_kgrid(self, dr, error, buffer, nufft_tol, reference)
      !> Grid instance
      type(moist_math_grid_3d_molecular_type), intent(inout) :: self
      !> Real-space spacing of the implied uniform box (bohr)
      real(wp), intent(in)    :: dr
      !> Error handling
      type(error_type), allocatable, intent(out)   :: error
      !> Optional bounding-box margin per side (bohr)
      real(wp), optional, intent(in)    :: buffer
      !> Optional FINUFFT tolerance; dominant cost control
      real(wp), optional, intent(in)    :: nufft_tol
      !> Optional center of symmetric period
      real(wp), optional, intent(in)    :: reference(3)

      real(wp) :: lo(3), hi(3), extent(3), box, marg
      integer  :: nk(3), d

      if (.not. allocated(self%xyz) .or. self%ngrid < 1) then
         call fatal_error(error, "molecular grid: build the grid before set_kgrid")
         return
      end if
      if (dr <= 0.0_wp) then
         call fatal_error(error, "molecular grid: set_kgrid needs dr > 0")
         return
      end if
      marg = default_kgrid_buffer
      if (present(buffer)) marg = buffer
      if (present(nufft_tol)) then
         if (nufft_tol <= 0.0_wp) then
            call fatal_error(error, "molecular grid: set_kgrid needs nufft_tol > 0")
            return
         end if
         self%nufft_tol = nufft_tol
      end if

      lo = minval(self%xyz, dim=2)
      hi = maxval(self%xyz, dim=2)
      extent = hi - lo
      if (present(reference)) extent = 2.0_wp*max(abs(hi - reference), abs(lo - reference))

      do d = 1, 3
         box = extent(d) + 2.0_wp*marg
         ! Check before CEILING to avoid integer overflow
         if (box/real(max_kgrid_modes_per_axis, wp) > dr) then
            call fatal_error(error, "molecular grid: auto-sized reciprocal grid is too "// &
               & "large (point cloud spans too wide for this dr); clamp the radial "// &
               & "extent (rcut_upper) when building the grid or use a coarser dr")
            return
         end if
         ! Smallest even N with N*dr covering extent
         nk(d) = ceiling(box/dr)
         nk(d) = max(nk(d), 2)
         if (mod(nk(d), 2) /= 0) nk(d) = nk(d) + 1
      end do

      if (any(nk > max_kgrid_modes_per_axis)) then
         call fatal_error(error, "molecular grid: auto-sized reciprocal grid is too "// &
            & "large (point cloud spans too wide for this dr); clamp the radial "// &
            & "extent (rcut_upper) when building the grid or use a coarser dr")
         return
      end if

      self%nkx = nk(1); self%nky = nk(2); self%nkz = nk(3)
      self%dkx = two_pi/(real(nk(1), wp)*dr)
      self%dky = two_pi/(real(nk(2), wp)*dr)
      self%dkz = two_pi/(real(nk(3), wp)*dr)
      self%npts_k = nk(1)*nk(2)*nk(3)
      ! Default phase reference `point(1)` matches Cartesian DFT
      self%centered_period = present(reference)
      self%kref_offset = 0.0_wp
      self%kref = self%xyz(:, 1)
      if (present(reference)) then
         self%kref = reference
         if (allocated(self%molecule)) then
            self%kref_offset = reference - sum(self%molecule%xyz, dim=2)/real(self%molecule%nat, wp)
         end if
      end if
      self%dr = dr
      self%kbuffer = marg
      self%has_kgrid = .true.
      ! Reciprocal changes invalidate earlier transform plans
      self%geom_generation = self%geom_generation + 1
   end subroutine molecular_grid_set_kgrid

   !> Reciprocal-space coordinate (1/bohr) of k-point j (1..npts_k)
   !>
   !> Invert CMCL index `j-1 = ax + nkx*(ay + nky*az)`
   !> Signed mode `m = a - N/2`; return `(m_x*dkx, m_y*dky, m_z*dkz)`
   !> Require `molecular_grid_set_kgrid` first
   !>
   !> @param[in]  self  Grid instance
   !> @param[in]  j     Flat k-point index (1..npts_k)
   pure function molecular_grid_kpoint(self, j) result(k)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> k-point index
      integer, intent(in) :: j
      !> Coordinate of k-point j
      real(wp) :: k(3)

      integer :: j0, ax, ay, az

      j0 = j - 1
      ax = mod(j0, self%nkx)
      ay = mod(j0/self%nkx, self%nky)
      az = j0/(self%nkx*self%nky)
      k = [real(ax - self%nkx/2, wp)*self%dkx, &
           real(ay - self%nky/2, wp)*self%dky, &
           real(az - self%nkz/2, wp)*self%dkz]
   end function molecular_grid_kpoint

   !> Allocate a transform bound to this grid
   !>
   !> @param[in]  self   Grid instance (must have the target attribute)
   !> @param[out] trafo  Allocated, bound transform engine
   !> @param[out] error  Set on allocation failure
   subroutine molecular_grid_new_trafo(self, trafo, error)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in), target :: self
      !> Allocated transform engine
      class(moist_math_grid_3d_trafo_type), allocatable, intent(out) :: trafo
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_3d_molecular_trafo_type), allocatable :: t
      integer :: stat

      allocate (t, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate molecular grid trafo")
         return
      end if
      t%grid => self
      t%nufft_tol = self%nufft_tol
      call move_alloc(t, trafo)
   end subroutine molecular_grid_new_trafo

   !> Bind an unprepared transform to `grid`
   !>
   !> Build k-grid, then call `trafo%prepare(nv)` for FINUFFT plans
   !> Grid must outlive transform
   !>
   !> @param[out] trafo      Initialized (unprepared) transform engine
   !> @param[in]  grid       Grid to transform on (must have the target attribute)
   !> @param[in]  nufft_tol  Optional requested FINUFFT relative tolerance; when
   !>                        absent the grid's own `nufft_tol` is used
   !> @param[in]  use_type12 Optional transform-family override: `.true.` (the
   !>                        default) uses the type-1/type-2 lattice path,
   !>                        `.false.` forces the reference type-3 path
   subroutine new_molecular_grid_trafo(trafo, grid, nufft_tol, use_type12)
      !> Initialized transform engine
      type(moist_math_grid_3d_molecular_trafo_type), intent(out) :: trafo
      !> Grid to transform on
      type(moist_math_grid_3d_molecular_type), intent(in), target :: grid
      !> Optional requested FINUFFT relative tolerance
      real(wp), optional, intent(in) :: nufft_tol
      !> Optional transform-family override
      logical, optional, intent(in) :: use_type12

      trafo%grid => grid
      trafo%nufft_tol = grid%nufft_tol
      if (present(nufft_tol)) trafo%nufft_tol = nufft_tol
      if (present(use_type12)) trafo%want_type12 = use_type12
   end subroutine new_molecular_grid_trafo

   !> Build FINUFFT plans for `ntrans = nv` transforms
   !>
   !> Both routes use `x_j = xyz(:,j) - kref` for Cartesian phase convention
   !>
   !> - Type 1/2 default: uniform-box dual lattice `k = m*dk`,
   !>   `dk = 2*pi/(nk*dr)`, `m = -N/2 .. N/2-1`
   !> - Phase identity `k.x = m*(dk*x)`; integer-mode periodicity
   !>   makes type 1/2 sums equal type 3 sums
   !> - Type 1 forward, type 2 backward; scaled coordinates `dk*x`
   !> - Type 1/2 fine grid depends on mode count; cheaper than type 3
   !> - CMCL `modeord = 0`, x fastest, preserves `f_k` ordering
   !> - Type 3 reference: explicit k coordinates, forward `iflag = -1`,
   !>   backward `iflag = +1`; supports elementwise path comparison
   !> - Fall back to type 3 outside portable `[-3pi, 3pi]` window
   !>
   !> Require configured k-grid; reprepare replaces existing plans
   !>
   !> @param[in,out] self    Trafo instance (bound to a grid)
   !> @param[in]     ntrans  Number of simultaneous transforms (= nv)
   !> @param[out]    error   Propagated error (unbound, no k-grid, FINUFFT, ...)
   subroutine molecular_trafo_prepare(self, ntrans, error)
      !> Trafo instance
      class(moist_math_grid_3d_molecular_trafo_type), intent(inout) :: self
      !> Number of simultaneous transforms
      integer, intent(in)    :: ntrans
      !> Error handling
      type(error_type), allocatable, intent(out)   :: error

      integer :: ax, ay, az, j, ngrid, npts_k, ier
      integer(int64) :: n_modes(3), m_npts, nk_npts, no_targets
      type(finufft_opts) :: opts
      real(c_double) :: dummy(1)
      real(wp) :: coord_span
      external :: finufft_makeplan, finufft_setpts, finufft_default_opts

      if (.not. associated(self%grid)) then
         call fatal_error(error, "molecular trafo: not bound to a grid")
         return
      end if
      if (.not. self%grid%has_kgrid) then
         call fatal_error(error, "molecular trafo: call molecular_grid_set_kgrid first")
         return
      end if
      if (ntrans < 1) then
         call fatal_error(error, "molecular trafo: ntrans must be >= 1")
         return
      end if

      ! Drop plans and scratch from previous preparation
      call molecular_trafo_free_plans(self)

      ngrid = self%grid%ngrid

      ! Shifted real-space coordinates
      ! `kref` shift preserves Cartesian FFT phase on both routes
      npts_k = self%grid%npts_k
      allocate (self%xj(ngrid), self%yj(ngrid), self%zj(ngrid))
      do j = 1, ngrid
         self%xj(j) = self%grid%xyz(1, j) - self%grid%kref(1)
         self%yj(j) = self%grid%xyz(2, j) - self%grid%kref(2)
         self%zj(j) = self%grid%xyz(3, j) - self%grid%kref(3)
      end do

      ! Decide the transform family
      ! Type 1/2 scales coordinates by per-axis dk: k.x = m.(dk*x)
      ! Spread <= 2*pi for a grid built by set_kgrid
      self%is_type12 = self%want_type12
      if (self%is_type12) then
         coord_span = max(maxval(abs(self%xj))*self%grid%dkx, &
                          maxval(abs(self%yj))*self%grid%dky, &
                          maxval(abs(self%zj))*self%grid%dkz)
         if (coord_span > type12_coord_limit) self%is_type12 = .false.
      end if

      if (self%is_type12) then
         do j = 1, ngrid
            self%xj(j) = self%xj(j)*self%grid%dkx
            self%yj(j) = self%yj(j)*self%grid%dky
            self%zj(j) = self%zj(j)*self%grid%dkz
         end do
      else
         ! Explicit nonuniform k coordinates, only needed by type 3
         allocate (self%xk(npts_k), self%yk(npts_k), self%zk(npts_k))
         j = 0
         do az = 0, self%grid%nkz - 1
            do ay = 0, self%grid%nky - 1
               do ax = 0, self%grid%nkx - 1
                  j = j + 1
                  self%xk(j) = real(ax - self%grid%nkx/2, c_double)*self%grid%dkx
                  self%yk(j) = real(ay - self%grid%nky/2, c_double)*self%grid%dky
                  self%zk(j) = real(az - self%grid%nkz/2, c_double)*self%grid%dkz
               end do
            end do
         end do
      end if
      allocate (self%cj(ngrid, ntrans))

      call finufft_default_opts(opts)
      ! CMCL modes: -N/2 .. N/2-1, x fastest
      opts%modeord = 0
      ! Explicit worker count: FINUFFT would otherwise take omp_get_max_threads
      opts%nthreads = int(self%grid%team_size(), c_int)

      m_npts = int(ngrid, int64)
      nk_npts = int(npts_k, int64)
      dummy = 0.0_c_double
      no_targets = 0_int64

      if (self%is_type12) then
         n_modes = [int(self%grid%nkx, int64), int(self%grid%nky, int64), &
                    int(self%grid%nkz, int64)]

         ! Forward: type 1, iflag = -1  =>  F(m) = sum_j c_j exp(-i m.(dk*x_j))
         call finufft_makeplan(1, 3, n_modes, -1, ntrans, real(self%nufft_tol, c_double), &
                               self%plan_r2k, opts, ier)
         if (ier /= 0) then
            call fatal_error(error, "molecular trafo: FINUFFT makeplan (forward) failed")
            call molecular_trafo_free_plans(self)
            return
         end if
         call finufft_setpts(self%plan_r2k, m_npts, self%xj, self%yj, self%zj, &
                             no_targets, dummy, dummy, dummy, ier)
         if (ier /= 0) then
            call fatal_error(error, "molecular trafo: FINUFFT setpts (forward) failed")
            call molecular_trafo_free_plans(self)
            return
         end if

         ! Backward: type 2, iflag = +1  =>  c_j = sum_m F(m) exp(+i m.(dk*x_j))
         call finufft_makeplan(2, 3, n_modes, 1, ntrans, real(self%nufft_tol, c_double), &
                               self%plan_k2r, opts, ier)
         if (ier /= 0) then
            call fatal_error(error, "molecular trafo: FINUFFT makeplan (backward) failed")
            call molecular_trafo_free_plans(self)
            return
         end if
         call finufft_setpts(self%plan_k2r, m_npts, self%xj, self%yj, self%zj, &
                             no_targets, dummy, dummy, dummy, ier)
         if (ier /= 0) then
            call fatal_error(error, "molecular trafo: FINUFFT setpts (backward) failed")
            call molecular_trafo_free_plans(self)
            return
         end if
      else
         n_modes = 0_int64  ! unused by FINUFFT type-3 plans

         ! Forward: type 3, iflag = -1  =>  F(k) = sum_j c_j exp(-i k.(r_j-kref))
         call finufft_makeplan(3, 3, n_modes, -1, ntrans, real(self%nufft_tol, c_double), &
                               self%plan_r2k, opts, ier)
         if (ier /= 0) then
            call fatal_error(error, "molecular trafo: FINUFFT makeplan (forward) failed")
            call molecular_trafo_free_plans(self)
            return
         end if
         call finufft_setpts(self%plan_r2k, m_npts, self%xj, self%yj, self%zj, &
                             nk_npts, self%xk, self%yk, self%zk, ier)
         if (ier /= 0) then
            call fatal_error(error, "molecular trafo: FINUFFT setpts (forward) failed")
            call molecular_trafo_free_plans(self)
            return
         end if

         ! Backward: type 3, iflag = +1  =>  c_j = sum_k F(k) exp(+i k.(r_j-kref))
         call finufft_makeplan(3, 3, n_modes, 1, ntrans, real(self%nufft_tol, c_double), &
                               self%plan_k2r, opts, ier)
         if (ier /= 0) then
            call fatal_error(error, "molecular trafo: FINUFFT makeplan (backward) failed")
            call molecular_trafo_free_plans(self)
            return
         end if
         call finufft_setpts(self%plan_k2r, nk_npts, self%xk, self%yk, self%zk, &
                             m_npts, self%xj, self%yj, self%zj, ier)
         if (ier /= 0) then
            call fatal_error(error, "molecular trafo: FINUFFT setpts (backward) failed")
            call molecular_trafo_free_plans(self)
            return
         end if
      end if

      ! Set readiness fields only after both plans succeed
      self%ntrans = ntrans
      self%geom_generation = self%grid%geom_generation
   end subroutine molecular_trafo_prepare

   !> Free FINUFFT plans and scratch
   !>
   !> Shared by reprepare and destroy; idempotent
   !>
   !> @param[in,out] self  Trafo instance
   subroutine molecular_trafo_free_plans(self)
      !> Trafo instance
      class(moist_math_grid_3d_molecular_trafo_type), intent(inout) :: self

      integer :: ier
      external :: finufft_destroy

      if (self%plan_r2k /= 0_int64) then
         call finufft_destroy(self%plan_r2k, ier)
         self%plan_r2k = 0_int64
      end if
      if (self%plan_k2r /= 0_int64) then
         call finufft_destroy(self%plan_k2r, ier)
         self%plan_k2r = 0_int64
      end if
      if (allocated(self%xj)) deallocate (self%xj)
      if (allocated(self%yj)) deallocate (self%yj)
      if (allocated(self%zj)) deallocate (self%zj)
      if (allocated(self%xk)) deallocate (self%xk)
      if (allocated(self%yk)) deallocate (self%yk)
      if (allocated(self%zk)) deallocate (self%zk)
      if (allocated(self%cj)) deallocate (self%cj)
      self%ntrans = 0
      self%is_type12 = .false.
   end subroutine molecular_trafo_free_plans

   !> Free plans and scratch; detach grid
   !>
   !> @param[in,out] self  Trafo instance
   subroutine molecular_trafo_destroy(self)
      !> Trafo instance
      class(moist_math_grid_3d_molecular_trafo_type), intent(inout) :: self

      call molecular_trafo_free_plans(self)
      self%grid => null()
   end subroutine molecular_trafo_destroy

   !> Forward NUFFT of all `nv` sites, real-space to reciprocal-space
   !>
   !> `f_r(ngrid, nv)` -> `f_k(npts_k, nv)` with strengths `w_j*f(r_j)`
   !> Matches Cartesian `dV*FFT`; requires `prepare(nv)` and matching `ntrans`
   !> Report precondition and FINUFFT failures through `error`
   !>
   !> @param[in,out] self   Trafo instance
   !> @param[in,out] f_r    Real-space field block, shape (ngrid, nv)
   !> @param[out]    f_k    Reciprocal-space field block, shape (npts_k, nv)
   !> @param[out]    error  Error handling
   subroutine molecular_trafo_fft_r2k(self, f_r, f_k, error)
      !> Trafo instance
      class(moist_math_grid_3d_molecular_trafo_type), intent(inout) :: self
      !> Real-space field block, shape (ngrid, nv)
      real(wp), intent(inout), contiguous, target :: f_r(:, :)
      !> Reciprocal-space field block, shape (npts_k, nv)
      complex(wp), intent(out), contiguous, target :: f_k(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: it, nv, ier
      external :: finufft_execute

      nv = size(f_r, 2)
      call molecular_trafo_check_ready(self, nv, size(f_r, 1), size(f_k, 1), error)
      if (allocated(error)) return
      ! Strengths c_j = w_j * f(r_j); real -> complex assignment zeroes imag
      do it = 1, nv
         self%cj(:, it) = self%grid%w*f_r(:, it)
      end do
      call finufft_execute(self%plan_r2k, self%cj, f_k, ier)
      if (ier /= 0) then
         call fatal_error(error, "molecular trafo: FINUFFT execute (forward) failed")
         return
      end if
   end subroutine molecular_trafo_fft_r2k

   !> Backward NUFFT of all `nv` sites, reciprocal-space to real-space
   !>
   !> `f_k(npts_k, nv)` -> `f_r(ngrid, nv)`
   !> Scale by `dkx*dky*dkz/(2*pi)^3`, equivalent to Cartesian `1/Vbox`
   !> Require `prepare(nv)`; report precondition and FINUFFT failures
   !>
   !> @param[in,out] self   Trafo instance
   !> @param[in,out] f_k    Reciprocal-space field block, shape (npts_k, nv)
   !> @param[out]    f_r    Real-space field block, shape (ngrid, nv)
   !> @param[out]    error  Error handling
   subroutine molecular_trafo_fft_k2r(self, f_k, f_r, error)
      !> Trafo instance
      class(moist_math_grid_3d_molecular_trafo_type), intent(inout) :: self
      !> Reciprocal-space field block, shape (npts_k, nv)
      complex(wp), intent(inout), contiguous, target :: f_k(:, :)
      !> Real-space field block, shape (ngrid, nv)
      real(wp), intent(out), contiguous, target :: f_r(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: it, nv, ier
      real(wp) :: scale
      external :: finufft_execute

      nv = size(f_k, 2)
      call molecular_trafo_check_ready(self, nv, size(f_r, 1), size(f_k, 1), error)
      if (allocated(error)) return
      ! FINUFFT `execute(plan, cj, fk)` argument roles depend on type
      ! Type 2: cj output, fk mode input; type 3: cj source, fk target
      ! Both paths store real-space values in `self%cj`
      if (self%is_type12) then
         call finufft_execute(self%plan_k2r, self%cj, f_k, ier)
      else
         call finufft_execute(self%plan_k2r, f_k, self%cj, ier)
      end if
      if (ier /= 0) then
         call fatal_error(error, "molecular trafo: FINUFFT execute (backward) failed")
         return
      end if
      scale = (self%grid%dkx*self%grid%dky*self%grid%dkz)/two_pi_cubed
      do it = 1, nv
         f_r(:, it) = scale*real(self%cj(:, it), wp)
      end do
   end subroutine molecular_trafo_fft_k2r

   !> Guard: trafo prepared, bound to its geometry, block sizes match
   !>
   !> Require prepared plans, matching geometry, and matching block sizes
   !>
   !> @param[in]  self   Trafo instance
   !> @param[in]  nv     Number of columns presented to the transform
   !> @param[in]  n_r    Leading extent of the caller's real-space block
   !> @param[in]  n_k    Leading extent of the caller's reciprocal-space block
   !> @param[out] error  Set when the trafo cannot transform this batch
   subroutine molecular_trafo_check_ready(self, nv, n_r, n_k, error)
      !> Trafo instance
      class(moist_math_grid_3d_molecular_trafo_type), intent(in) :: self
      !> Batch width of this call
      integer, intent(in) :: nv
      !> Leading extent of the caller's real-space block
      integer, intent(in) :: n_r
      !> Leading extent of the caller's reciprocal-space block
      integer, intent(in) :: n_k
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (self%ntrans < 1 .or. .not. allocated(self%cj)) then
         call fatal_error(error, "molecular trafo: call prepare(nv) before transforming")
         return
      end if
      if (self%geom_generation /= self%grid%geom_generation) then
         call fatal_error(error, &
            & "molecular trafo: grid geometry changed since prepare; destroy and recreate the trafo")
         return
      end if
      if (nv /= self%ntrans) then
         call fatal_error(error, "molecular trafo: batch width does not match prepared ntrans")
         return
      end if
      if (n_r /= self%grid%ngrid .or. n_k /= self%grid%npts_k .or. size(self%cj, 1) /= self%grid%ngrid) then
         call fatal_error(error, "molecular trafo: field block size does not match the prepared grid")
         return
      end if
   end subroutine molecular_trafo_check_ready

   !> Write a short grid summary to the given unit
   !>
   !> @param[in] self  Grid instance
   !> @param[in] unit  Optional output unit (default `output_unit`)
   subroutine molecular_grid_info(self, unit)
      !> Grid instance
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Output unit (default `output_unit`)
      integer, optional, intent(in) :: unit

      integer :: iunit, iat, nat

      iunit = output_unit
      if (present(unit)) iunit = unit

      write (iunit, "(a)") "moist moist_math_grid_3d_molecular_type"
      write (iunit, "(a,i0)") "  total points   : ", self%ngrid
      if (.not. allocated(self%xyz)) then
         write (iunit, "(a)") "  (grid is uninitialised)"
         return
      end if
      write (iunit, "(a,es12.4,a,es12.4)") &
         & "  weight range   : min = ", minval(self%w), &
         & ", max = ", maxval(self%w)
      write (iunit, "(a,es14.6)") "  sum(weights)   : ", sum(self%w)
      nat = size(self%nrad_per_atom)
      write (iunit, "(a,i0)") "  atoms          : ", nat
      write (iunit, "(a)") "  per-atom counts (atom, ngrid, nrad, nang):"
      do iat = 1, nat
         write (iunit, "(4x,i6,3x,i8,3x,i5,3x,i5)") iat, &
            & self%atom_offset(iat + 1) - self%atom_offset(iat), &
            & self%nrad_per_atom(iat), self%nang_per_atom(iat)
      end do
   end subroutine molecular_grid_info

end module moist_math_grid_3d_molecular
