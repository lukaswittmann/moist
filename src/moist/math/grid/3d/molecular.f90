!> Atom-centered molecular integration grid
module moist_math_grid_3d_molecular
   use, intrinsic :: iso_fortran_env, only: output_unit, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use moist_data_atomicrad, only: covalent_rad
   use moist_math_quadrature_lebedev, only: get_angular_grid, grid_size, lebedev_order_from_num
   use moist_math_quadrature_chebyshev, only: chebyshev2_radii
   use moist_math_quadrature_handymod, only: handymod_radii, handymod_midpoint_radii
   use moist_math_quadrature_becke, only: becke_weights
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type, &
                                      integrand_3d
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
   use finufft_mod, only: finufft_opts
!$ use omp_lib, only: omp_get_max_threads, omp_in_parallel
   implicit none(type, external)
   private

   public :: moist_math_grid_3d_molecular_type
   public :: new_molecular_grid
   public :: new_molecular_grid_uniform
   public :: new_molecular_grid_uniform_handymod
   public :: new_molecular_grid_uniform_qc_handymod
   public :: lebedev_degree_from_num
   public :: default_grid_sizes
   public :: molecular_grid_set_kgrid
   public :: moist_math_grid_3d_molecular_trafo_type, new_molecular_grid_trafo
   public :: integrand_3d
   public :: default_nufft_tol
   public :: lebedev_order_from_arc

   !> Default requested FINUFFT relative tolerance
   real(wp), parameter :: default_nufft_tol = 1.0e-10_wp

   !> Lebedev sizes carrying negative quadrature weights
   !>
   !> The per-shell angular rule skips them: a negative weight makes the
   !> shell integral non-monotone under refinement, which is exactly the
   !> property the rule relies on when it trades resolution for cost
   integer, parameter :: lebedev_negative_weight_sizes(3) = [74, 230, 266]

   !> Default floor on the per-shell Lebedev point count
   !>
   !> 110 is the smallest order any converging molecular 3D-RISM
   !> configuration has ever used; the inner shells it over-resolves cost
   !> well under 1% of the grid, so the floor is kept in tested territory
   !> rather than pushed to the arc target
   integer, parameter :: default_shell_nang_min = 110

   !> Default cap on the per-shell Lebedev point count
   !>
   !> The largest grid the Lebedev-Laikov tables provide
   integer, parameter :: default_shell_nang_max = 5810

   !> Default weight-pruning threshold, bohr^3
   !>
   !> Points with |w| below this are dropped from the final grid
   real(wp), parameter :: default_wthr = 1.0e-14_wp

   !> Algebraic Lebedev degrees corresponding to the supported angular grids
   integer, parameter :: lebedev_degree_table(32) = [ &
      &   3, 5, 7, 9, 11, 13, 15, 17, &
      &  19, 21, 23, 25, 27, 29, 31, 35, &
      &  41, 47, 53, 59, 65, 71, 77, 83, &
      &  89, 95, 101, 107, 113, 119, 125, 131]

   !> Two pi, used for the reciprocal-grid spacing dk = 2*pi/(N*dr)
   real(wp), parameter :: two_pi = 8.0_wp*atan(1.0_wp)
   !> (2*pi)^3, the inverse-Fourier-transform normalisation denominator
   real(wp), parameter :: two_pi_cubed = two_pi*two_pi*two_pi
   !> Largest magnitude of a scaled type-1/type-2 coordinate handed to FINUFFT
   !>
   !> The bundled FINUFFT (2.6) folds arbitrary coordinates, but older
   !> releases restrict them to `[-3pi, 3pi]`; staying inside that window
   !> keeps the fast path portable, and anything outside it silently falls
   !> back to the type-3 path instead of risking a wrong answer
   real(wp), parameter :: type12_coord_limit = 3.0_wp*(0.5_wp*two_pi)
   !> Default extra real-space margin, bohr, on each side of the point box
   !>
   !> Added when auto-sizing the reciprocal grid; keeps the shifted
   !> coordinates fed to FINUFFT clear of its periodic wrap
   real(wp), parameter :: default_kgrid_buffer = 2.0_wp
   !> Safety cap on the auto-sized reciprocal-grid modes per axis
   !>
   !> A molecular grid keeps far, low-weight radial shells, so without a
   !> radial clamp the bounding box (and hence N = box/dr) can explode and
   !> exhaust memory in the NUFFT
   !> Exceeding this raises a clear error instead of being OOM-killed
   integer, parameter :: max_kgrid_modes_per_axis = 256

   !> Molecular construction recipes
   integer, parameter :: recipe_element_defaults = 1, recipe_uniform = 2
   !> Integer and midpoint HandyMod construction recipes
   integer, parameter :: recipe_handymod = 3, recipe_qc_handymod = 4

   !> Atom-centered molecular integration grid
   !>
   !> - Chebyshev-2 radial x Lebedev angular x Becke fuzzy-cell partitioning
   !> - Coordinates, count and volume fields inherited from the domain bases
   !> - No transform handles; destroy borrowed engines before geometry
   !>   mutation
   type, extends(moist_math_grid_3d_type) :: moist_math_grid_3d_molecular_type
      !> Radial midpoint count per atom
      integer :: nrad = 50
      !> Algebraic Lebedev degree
      integer :: lebedev_degree = 17
      !> Inner radial boundary, bohr
      real(wp) :: rmin = 0.0_wp
      !> Outer radial boundary, bohr
      real(wp) :: rmax = 10.0_wp
      !> HandyMod mapping parameter
      real(wp) :: m = 2.0_wp
      !> Becke polynomial iteration count
      integer :: becke_k = 3
      !> Optional SSF cell parameter, replacing the iterated polynomial
      real(wp), allocatable :: ssf_a
      !> Bare Becke weight pruning threshold
      real(wp) :: pruning_threshold = 0.0_wp
      !> Real-space spacing used to size the initial reciprocal grid, bohr
      real(wp) :: dr = 0.5_wp
      !> Period margin on each side of the point cloud, bohr
      real(wp) :: kbuffer = 2.0_wp
      !> Publish weight-dependent Gaussian widths instead of point potentials
      logical :: gaussian = .false.
      !> Width scale, xi0 = xi0_factor / w**(1/3), for positive weights
      real(wp) :: xi0_factor = 1.0_wp
      !> Angular count for Chebyshev and integer HandyMod construction
      integer :: nang = 26
      !> Optional angular band edges
      real(wp), allocatable :: arc_r(:)
      !> Optional angular band spacings
      real(wp), allocatable :: arc_a(:)
      !> Angular count floor for band-based construction
      integer :: nang_min = default_shell_nang_min
      !> Angular count cap for band-based construction
      integer :: nang_max = default_shell_nang_max
      !> Radial/angular construction selected by the constructor
      integer, private :: radial_rule = recipe_qc_handymod
      !> Configure reciprocal geometry on the first domain update
      logical, private :: auto_kgrid = .true.
      !> CSR-style atom ownership, shape (nat+1)
      !>
      !> Points owned by atom i stored contiguously at indices
      !> [atom_offset(i), atom_offset(i+1)-1]; atom_offset(1) = 1,
      !> atom_offset(nat+1) = ngrid+1
      integer, allocatable :: atom_offset(:)
      !> Per-atom radial size actually used, shape (nat)
      integer, allocatable :: nrad_per_atom(:)
      !> Per-atom Lebedev point count actually used, shape (nat)
      integer, allocatable :: nang_per_atom(:)

      ! --- Uniform reciprocal (k) grid for the NUFFT ---
      !> Mode counts per axis
      !>
      !> Built by `molecular_grid_set_kgrid`; until then `has_kgrid` is
      !> false, `npts_k` is 0, and no transform can be issued
      !> K-grid is a uniform box of `nkx*nky*nkz` modes with per-axis
      !> spacing `dk*`, enclosing all real-space points
      !> Mode ordering follows FINUFFT's default CMCL layout
      !> (`modeord = 0`): kx fastest, frequencies -N/2 .. N/2-1
      integer :: nkx = 0, nky = 0, nkz = 0
      !> Reciprocal-space spacings (1/bohr)
      real(wp) :: dkx = 0.0_wp, dky = 0.0_wp, dkz = 0.0_wp
      !> Phase reference: real-space point the transform phases are measured from
      !>
      !> Direct grids default to `point(1)`; default domain updates select
      !> the solute centroid
      real(wp) :: kref(3) = 0.0_wp
      !> Requested FINUFFT relative tolerance for transforms on this grid
      !>
      !> Set through `molecular_grid_set_kgrid(..., nufft_tol=...)`; the
      !> trafo picks it up in `prepare` unless overridden per trafo
      !> See `default_nufft_tol` for the cost/accuracy trade this controls
      real(wp) :: nufft_tol = default_nufft_tol
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
      !> Bumped whenever the committed point cloud or reciprocal grid
      !> changes (`move_molecular_geometry`, `molecular_grid_set_kgrid`,
      !> `destroy`)
      !> A trafo captures this at `prepare` and `check_ready` rejects a
      !> mismatch, so a trafo cannot silently run its plans against a
      !> geometry they were not built for
      integer, private :: geom_generation = 0
   contains
      procedure :: validate => validate_molecular_grid
      procedure :: update => update_molecular
      procedure :: rebuild => rebuild_molecular
      procedure :: kind_name => molecular_kind_name
      procedure :: has_geometry_dependent_xi0 => molecular_xi0_dependent
      !> Reciprocal-space coordinate of k-point j (1..npts_k; needs k-grid set)
      procedure :: kpoint => molecular_grid_kpoint
      !> Free all array components (idempotent)
      procedure :: destroy => molecular_grid_destroy
      !> Print a short summary to the given unit (default output_unit)
      procedure :: info => molecular_grid_info
      !> Allocate the matching transform engine bound to this grid
      procedure :: new_trafo => molecular_grid_new_trafo
   end type moist_math_grid_3d_molecular_type

   !> Fourier transform engine for a molecular grid
   !>
   !> A molecular grid is non-uniform, so its transform is a non-uniform FFT
   !> (NUFFT), backed by FINUFFT.  Forward
   !> (real -> reciprocal) carries strengths `w_j * f(r_j)` (the quadrature
   !> weights baked in, mirroring the Cartesian `dV * FFT` convention);
   !> backward is scaled by `dkx*dky*dkz/(2*pi)^3`.  Both use the grid's phase
   !> reference `kref` so the result matches the Cartesian/`potential_3d`
   !> convention exactly
   !>
   !> Because the k-points are the FFT dual lattice of a uniform box (see
   !> `molecular_trafo_prepare`), the pair is issued as FINUFFT type 1
   !> (forward) and type 2 (backward) rather than the general nonuniform ->
   !> nonuniform type 3; `want_type12 = .false.` selects the type-3 reference
   !> route, which computes the same sums at higher cost
   !>
   !> `prepare(ntrans)` builds two guru plans sized for `ntrans = nv`
   !> simultaneous transforms; the molecular points are sorted once at that
   !> point and reused for every later `fft_r2k`/`fft_k2r`.  The plans hold
   !> mutable internal scratch, so a prepared trafo must be executed from a
   !> single thread (the batched call threads `nv` internally) rather than
   !> shared across OpenMP threads
   !> The plans also fix the worker count at `prepare` (one inside an OpenMP
   !> region, otherwise `omp_get_max_threads()`); a later thread-count change
   !> takes effect only after the next `prepare`
   type, extends(moist_math_grid_3d_trafo_type) :: moist_math_grid_3d_molecular_trafo_type
      !> Grid this trafo transforms on (not owned; must outlive the trafo)
      class(moist_math_grid_3d_molecular_type), pointer :: grid => null()
      !> Number of simultaneous transforms baked into the plans (0 = unprepared)
      integer :: ntrans = 0
      !> Grid's `geom_generation` at the end of a successful `prepare`
      !>
      !> `check_ready` compares it against the grid's current value to catch
      !> a trafo whose plans outlived a geometry change
      integer :: geom_generation = -1
      !> Requested FINUFFT relative tolerance actually used by the plans
      !>
      !> Seeded from the grid's `nufft_tol` unless `new_molecular_grid_trafo`
      !> was given an explicit override
      real(wp) :: nufft_tol = default_nufft_tol
      !> Request the type-1/type-2 fast path (the default)
      !>
      !> `.false.` forces the reference type-3 path; the two are
      !> mathematically identical (see `molecular_trafo_prepare`) and exist
      !> side by side so the fast path can be A/B-checked elementwise inside
      !> one process
      logical :: want_type12 = .true.
      !> Whether the prepared plans are actually type 1/2 (`.false.` = type 3)
      !>
      !> `prepare` falls back to type 3 if the scaled coordinates would leave
      !> FINUFFT's portable `[-3pi, 3pi]` input window
      logical :: is_type12 = .false.
      !> Opaque FINUFFT plan handles, real-space to reciprocal-space (0 = none)
      integer(int64) :: plan_r2k = 0_int64
      !> Opaque FINUFFT plan handle, reciprocal-space to real-space (0 = none)
      integer(int64) :: plan_k2r = 0_int64
      !> Shifted molecular-grid coordinates, length ngrid
      !>
      !> Kept alive for the plan lifetime
      !> Type 3: `xj(j) = xyz(1,j) - kref(1)` in bohr
      !> Type 1/2: the same shifted coordinate scaled by `dkx` into
      !> FINUFFT's radian convention, `xj(j) = (xyz(1,j) - kref(1))*dkx`
      real(c_double), allocatable :: xj(:), yj(:), zj(:)
      !> Reciprocal-grid target/source coordinates, length npts_k
      !>
      !> Allocated only on the type-3 path; type 1/2 addresses the lattice
      !> by its mode indices and needs no explicit k coordinates
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

   !> Validate settings before constructing quadrature or reciprocal modes
   !>
   !> @param[in] self Grid construction settings
   !> @param[out] error Invalid setting
   subroutine validate_molecular_grid(self, error)
      !> Grid construction settings
      class(moist_math_grid_3d_molecular_type), intent(in) :: self
      !> Invalid setting
      type(error_type), allocatable, intent(out) :: error

      if (.not. all(ieee_is_finite([self%rmin, self%rmax, self%m, self%pruning_threshold, &
                                    self%dr, self%kbuffer, self%nufft_tol, self%xi0_factor]))) then
         call fatal_error(error, "molecular domain: settings must be finite")
      else if (self%nrad < 1 .or. self%lebedev_degree < 1 .or. self%becke_k < 0) then
         call fatal_error(error, "molecular domain: invalid quadrature order or Becke iteration count")
      else if (self%rmin < 0.0_wp .or. self%rmax <= self%rmin .or. self%m <= 0.0_wp) then
         call fatal_error(error, "molecular domain: invalid HandyMod radii or mapping parameter")
      else if (self%radial_rule == recipe_qc_handymod .and. &
         & self%m >= log(self%rmax - self%rmin + 1.0_wp)/log(2.0_wp)) then
         call fatal_error(error, "molecular domain: HandyMod requires rmax - rmin > 2**m - 1")
      else if (self%pruning_threshold < 0.0_wp .or. self%pruning_threshold >= 1.0_wp) then
         call fatal_error(error, "molecular domain: pruning_threshold must be in [0, 1)")
      else if (self%dr <= 0.0_wp .or. self%kbuffer < 0.0_wp .or. self%nufft_tol <= 0.0_wp) then
         call fatal_error(error, "molecular domain: invalid dr, kbuffer or nufft_tol")
      else if (self%xi0_factor <= 0.0_wp) then
         call fatal_error(error, "molecular domain: xi0_factor must be positive")
      end if
      if (allocated(error)) return
      if (self%radial_rule == recipe_handymod) then
         if (self%m > real(huge(0), wp) .or. self%m /= aint(self%m)) then
            call fatal_error(error, "molecular domain: integer HandyMod requires a representable integer m")
            return
         end if
      end if
      if (allocated(self%ssf_a)) then
         if (.not. ieee_is_finite(self%ssf_a)) then
            call fatal_error(error, "molecular domain: ssf_a must be finite")
         else if (self%ssf_a <= 0.0_wp .or. self%ssf_a > 1.0_wp) then
            call fatal_error(error, "molecular domain: ssf_a must be in (0, 1]")
         end if
      end if
   end subroutine validate_molecular_grid

   !> Retain direct-constructor settings for subsequent domain updates
   !>
   !> @param[in,out] self Constructed grid
   !> @param[in] mol Molecular structure
   !> @param[in] flavor Constructor selection
   !> @param[in] becke_k Becke iteration count
   !> @param[in] becke_ssf_a Optional SSF cell parameter
   !> @param[in] becke_thr Becke pruning threshold
   !> @param[in] nrad Uniform radial size
   !> @param[in] nang Angular point count
   !> @param[in] lebedev_degree Algebraic angular degree
   !> @param[in] rmin Lower radial bound
   !> @param[in] rmax Upper radial bound
   !> @param[in] handymod_m Integer HandyMod exponent
   !> @param[in] m Real HandyMod parameter
   !> @param[in] arc_r Angular band edges
   !> @param[in] arc_a Angular band spacings
   !> @param[in] nang_min Angular count floor
   !> @param[in] nang_max Angular count cap
   subroutine remember_molecular_settings(self, mol, flavor, becke_k, becke_ssf_a, becke_thr, &
      & nrad, nang, rmin, rmax, handymod_m, m, arc_r, arc_a, nang_min, nang_max, lebedev_degree)
      !> Constructed grid
      type(moist_math_grid_3d_molecular_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Constructor selection
      integer, intent(in) :: flavor
      !> Becke iteration count
      integer, optional, intent(in) :: becke_k
      !> SSF cell parameter
      real(wp), optional, intent(in) :: becke_ssf_a
      !> Becke pruning threshold
      real(wp), optional, intent(in) :: becke_thr
      !> Uniform radial size
      integer, optional, intent(in) :: nrad
      !> Angular point count
      integer, optional, intent(in) :: nang
      !> Algebraic angular degree
      integer, optional, intent(in) :: lebedev_degree
      !> Lower radial bound
      real(wp), optional, intent(in) :: rmin
      !> Upper radial bound
      real(wp), optional, intent(in) :: rmax
      !> Integer HandyMod exponent
      integer, optional, intent(in) :: handymod_m
      !> Real HandyMod parameter
      real(wp), optional, intent(in) :: m
      !> Angular band edges
      real(wp), optional, intent(in) :: arc_r(:)
      !> Angular band spacings
      real(wp), optional, intent(in) :: arc_a(:)
      !> Angular count floor
      integer, optional, intent(in) :: nang_min
      !> Angular count cap
      integer, optional, intent(in) :: nang_max

      self%molecule = mol
      self%radial_rule = flavor
      self%auto_kgrid = .false.
      self%rmin = 0.0_wp
      self%rmax = huge(1.0_wp)
      if (present(becke_k)) self%becke_k = becke_k
      if (present(becke_ssf_a)) self%ssf_a = becke_ssf_a
      if (present(becke_thr)) self%pruning_threshold = becke_thr
      if (present(nrad)) self%nrad = nrad
      if (present(nang)) self%nang = nang
      if (present(lebedev_degree)) self%lebedev_degree = lebedev_degree
      if (present(rmin)) self%rmin = rmin
      if (present(rmax)) self%rmax = rmax
      if (present(handymod_m)) self%m = real(handymod_m, wp)
      if (present(m)) self%m = m
      if (present(arc_r)) self%arc_r = arc_r
      if (present(arc_a)) self%arc_a = arc_a
      if (present(nang_min)) self%nang_min = nang_min
      if (present(nang_max)) self%nang_max = nang_max
   end subroutine remember_molecular_settings

   !> Commit a successfully constructed molecular discretization
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
      ! The committed point cloud changed (a fresh candidate was built and
      ! moved in), so any trafo prepared against the old one is now stale
      dest%geom_generation = dest%geom_generation + 1
      call move_alloc(source%xyz, dest%xyz)
      call move_alloc(source%w, dest%w)
      call move_alloc(source%owner, dest%owner)
      call move_alloc(source%xi0, dest%xi0)
      call move_alloc(source%atom_offset, dest%atom_offset)
      call move_alloc(source%nrad_per_atom, dest%nrad_per_atom)
      call move_alloc(source%nang_per_atom, dest%nang_per_atom)
   end subroutine move_molecular_geometry

   !> Move atomic grids, recompute Becke weights and guard the fixed period
   !>
   !> @param[in,out] self Domain instance
   !> @param[in] mol Solute structure
   !> @param[out] error Invalid geometry or insufficient period
   subroutine update_molecular(self, mol, error)
      !> Domain instance
      class(moist_math_grid_3d_molecular_type), intent(inout) :: self
      !> Solute structure
      type(structure_type), intent(in) :: mol
      !> Invalid geometry or insufficient period
      type(error_type), allocatable, intent(out) :: error

      if (mol%nat < 1) then
         call fatal_error(error, "molecular domain: at least one solute atom is required")
         return
      end if
      if (.not. all(ieee_is_finite(mol%xyz))) then
         call fatal_error(error, "molecular domain: solute coordinates must be finite")
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
   !> Fixed orders regenerate the same local quadrature points
   !> Reapplying threshold pruning permits point-count changes at the
   !> documented cutoff
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
      !> Point owner index
      integer :: stat
      !> Whether to retain the previous reciprocal grid
      logical :: preserve_period
      !> Fourier period factor
      real(wp), parameter :: two_pi = 8.0_wp*atan(1.0_wp)

      call self%validate(error)
      if (allocated(error)) return
      allocate (grid, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "molecular domain: cannot allocate grid")
         return
      end if
      associate (s => self, mol => self%molecule)
         center = sum(mol%xyz, dim=2)/real(mol%nat, wp)
         select case (s%radial_rule)
         case (recipe_element_defaults)
            call new_molecular_grid(grid, mol, error, s%becke_k, s%ssf_a, s%pruning_threshold)
         case (recipe_uniform)
            call new_molecular_grid_uniform(grid, mol, s%nrad, s%nang, error, &
               & s%rmin, s%rmax, s%becke_k, s%ssf_a, s%pruning_threshold)
         case (recipe_handymod)
            call new_molecular_grid_uniform_handymod(grid, mol, s%nrad, s%nang, error, &
               & s%rmin, s%rmax, int(s%m), s%becke_k, s%ssf_a, s%pruning_threshold)
         case (recipe_qc_handymod)
            call new_molecular_grid_uniform_qc_handymod(grid, mol, s%nrad, s%lebedev_degree, error, &
               & s%rmin, s%rmax, s%m, s%becke_k, s%ssf_a, s%pruning_threshold, &
               & s%arc_r, s%arc_a, s%nang_min, s%nang_max)
         case default
            call fatal_error(error, "molecular domain: unknown quadrature recipe")
         end select
      end associate
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

   !> Literature-standard per-element default grid sizes (Pople "fine"/"ultrafine" style)
   !>
   !>   H, He            -> (50, 302)
   !>   Li..Ne (row 2)   -> (75, 302)
   !>   Na..Ar (row 3)   -> (75, 434)
   !>   K  and beyond    -> (99, 590)
   !>
   !> `nang` is the raw Lebedev point count; it is always one of the
   !> supported sizes
   !>
   !> @param[in]  iz    Atomic number
   !> @param[out] nrad  Number of radial points
   !> @param[out] nang  Raw Lebedev point count
   pure subroutine default_grid_sizes(iz, nrad, nang)
      !> Atomic number
      integer, intent(in)  :: iz
      !> Number of radial points
      integer, intent(out) :: nrad
      !> Raw Lebedev point count
      integer, intent(out) :: nang

      if (iz <= 2) then
         nrad = 50; nang = 302
      else if (iz <= 10) then
         nrad = 75; nang = 302
      else if (iz <= 18) then
         nrad = 75; nang = 434
      else
         nrad = 99; nang = 590
      end if
   end subroutine default_grid_sizes

   !> Build a molecular grid using element-dependent default sizes
   !>
   !> @param[out] self         Initialised molecular grid
   !> @param[in]  mol          Molecular structure (mctc_io `structure_type`)
   !> @param[out] error        Propagated error (invalid grid sizes, alloc, ...)
   !> @param[in]  becke_k      Optional Becke polynomial stiffness (default 3)
   !> @param[in]  becke_ssf_a  Optional SSF cell parameter, replaces the polynomial
   !> @param[in]  becke_thr    Optional bare-Becke-weight prune threshold (default 0)
   subroutine new_molecular_grid(self, mol, error, becke_k, becke_ssf_a, becke_thr)
      !> Grid to initialise
      type(moist_math_grid_3d_molecular_type), intent(out) :: self
      !> Molecular structure
      type(structure_type), intent(in)  :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional Becke cutoff-polynomial stiffness `k` (default 3)
      !>
      !> Shared across the grid constructors, so absent must mean exactly
      !> Becke's k = 3
      integer, optional, intent(in)  :: becke_k
      !> Optional SSF cell-function parameter; replaces the Becke polynomial
      real(wp), optional, intent(in)  :: becke_ssf_a
      !> Optional bare-Becke-weight prune threshold
      !>
      !> Points whose partition weight is below this are dropped as owned
      !> by another atom; default 0 prunes nothing beyond the total-weight
      !> threshold `default_wthr`
      real(wp), optional, intent(in)  :: becke_thr

      integer, allocatable :: nrad_atom(:), nang_atom(:)
      integer :: iat, iz

      allocate (nrad_atom(mol%nat), nang_atom(mol%nat))
      do iat = 1, mol%nat
         iz = mol%num(mol%id(iat))
         call default_grid_sizes(iz, nrad_atom(iat), nang_atom(iat))
      end do

      call build_molecular_grid(self, mol, nrad_atom, nang_atom, &
         & default_wthr, error, becke_k=becke_k, becke_ssf_a=becke_ssf_a, becke_thr=becke_thr)
      if (allocated(error)) return
      call remember_molecular_settings(self, mol, recipe_element_defaults, becke_k, becke_ssf_a, becke_thr)
   end subroutine new_molecular_grid

   !> Build a molecular grid with uniform `(nrad, nang)` per atom
   !>
   !> @param[out] self         Initialised molecular grid
   !> @param[in]  mol          Molecular structure
   !> @param[in]  nrad         Number of radial points per atom (>= 1)
   !> @param[in]  nang         Raw Lebedev point count per atom (must be one
   !>                          of the supported sizes: 6, 14, 26, ..., 5810)
   !> @param[out] error        Propagated error
   !> @param[in]  rmin         Optional minimum radial shell radius (bohr)
   !> @param[in]  rmax         Optional maximum radial shell radius (bohr)
   !> @param[in]  becke_k      Optional Becke polynomial stiffness (default 3)
   !> @param[in]  becke_ssf_a  Optional SSF cell parameter, replaces the polynomial
   !> @param[in]  becke_thr    Optional bare-Becke-weight prune threshold (default 0)
   subroutine new_molecular_grid_uniform(self, mol, nrad, nang, error, rmin, rmax, &
         & becke_k, becke_ssf_a, becke_thr)
      !> Grid to initialise
      type(moist_math_grid_3d_molecular_type), intent(out) :: self
      !> Molecular structure
      type(structure_type), intent(in)  :: mol
      !> Number of radial points per atom
      integer, intent(in)  :: nrad
      !> Raw Lebedev point count per atom
      integer, intent(in)  :: nang
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Minimum radius (bohr)
      real(wp), optional, intent(in)  :: rmin
      !> Maximum radius (bohr)
      real(wp), optional, intent(in)  :: rmax
      !> Optional Becke cutoff-polynomial stiffness `k` (default 3)
      !>
      !> Shared across the grid constructors, so absent must mean exactly
      !> Becke's k = 3
      integer, optional, intent(in)  :: becke_k
      !> Optional SSF cell-function parameter; replaces the Becke polynomial
      real(wp), optional, intent(in)  :: becke_ssf_a
      !> Optional bare-Becke-weight prune threshold
      !>
      !> Points whose partition weight is below this are dropped as owned
      !> by another atom; default 0 prunes nothing beyond the total-weight
      !> threshold `default_wthr`
      real(wp), optional, intent(in)  :: becke_thr

      integer, allocatable :: nrad_atom(:), nang_atom(:)

      if (nrad < 1) then
         call fatal_error(error, "molecular grid: nrad must be >= 1")
         return
      end if

      allocate (nrad_atom(mol%nat), nang_atom(mol%nat))
      nrad_atom = nrad
      nang_atom = nang

      call build_molecular_grid(self, mol, nrad_atom, nang_atom, &
         & default_wthr, error, rmin=rmin, rmax=rmax, becke_k=becke_k, becke_ssf_a=becke_ssf_a, becke_thr=becke_thr)
      if (allocated(error)) return
      call remember_molecular_settings(self, mol, recipe_uniform, becke_k, becke_ssf_a, becke_thr, &
         & nrad=nrad, nang=nang, rmin=rmin, rmax=rmax)
   end subroutine new_molecular_grid_uniform

   !> Build a molecular grid with uniform `(nrad, nang)` per atom and HandyMod radii
   !>
   !> @param[out] self         Initialised molecular grid
   !> @param[in]  mol          Molecular structure
   !> @param[in]  nrad         Number of radial points per atom (>= 1)
   !> @param[in]  nang         Raw Lebedev point count per atom (must be supported)
   !> @param[out] error        Propagated error
   !> @param[in]  rmin         Minimum HandyMod radius (bohr)
   !> @param[in]  rmax         Maximum HandyMod radius (bohr)
   !> @param[in]  m            HandyMod exponent (>= 1)
   !> @param[in]  becke_k      Optional Becke polynomial stiffness (default 3)
   !> @param[in]  becke_ssf_a  Optional SSF cell parameter, replaces the polynomial
   !> @param[in]  becke_thr    Optional bare-Becke-weight prune threshold (default 0)
   subroutine new_molecular_grid_uniform_handymod(self, mol, nrad, nang, error, rmin, rmax, m, &
         & becke_k, becke_ssf_a, becke_thr)
      !> Grid to initialise
      type(moist_math_grid_3d_molecular_type), intent(out) :: self
      !> Molecular structure
      type(structure_type), intent(in)  :: mol
      !> Number of radial points per atom
      integer, intent(in)  :: nrad
      !> Raw Lebedev point count per atom
      integer, intent(in)  :: nang
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Minimum HandyMod radius (bohr)
      real(wp), intent(in)  :: rmin
      !> Maximum HandyMod radius (bohr)
      real(wp), intent(in)  :: rmax
      !> HandyMod exponent
      integer, intent(in)  :: m
      !> Optional Becke cutoff-polynomial stiffness `k` (default 3)
      !>
      !> Shared across the grid constructors, so absent must mean exactly
      !> Becke's k = 3
      integer, optional, intent(in)  :: becke_k
      !> Optional SSF cell-function parameter; replaces the Becke polynomial
      real(wp), optional, intent(in)  :: becke_ssf_a
      !> Optional bare-Becke-weight prune threshold
      !>
      !> Points whose partition weight is below this are dropped as owned
      !> by another atom; default 0 prunes nothing beyond the total-weight
      !> threshold `default_wthr`
      real(wp), optional, intent(in)  :: becke_thr

      integer, allocatable :: nrad_atom(:), nang_atom(:)

      if (nrad < 1) then
         call fatal_error(error, "molecular grid: nrad must be >= 1")
         return
      end if
      if (rmax <= rmin) then
         call fatal_error(error, "molecular grid HandyMod: rmax must be greater than rmin")
         return
      end if
      if (m < 1) then
         call fatal_error(error, "molecular grid HandyMod: m must be >= 1")
         return
      end if

      allocate (nrad_atom(mol%nat), nang_atom(mol%nat))
      nrad_atom = nrad
      nang_atom = nang

      call build_molecular_grid(self, mol, nrad_atom, nang_atom, &
         & default_wthr, error, rmin=rmin, rmax=rmax, handymod_m=m, becke_k=becke_k, becke_ssf_a=becke_ssf_a, becke_thr=becke_thr)
      if (allocated(error)) return
      call remember_molecular_settings(self, mol, recipe_handymod, becke_k, becke_ssf_a, becke_thr, &
         & nrad=nrad, nang=nang, rmin=rmin, rmax=rmax, handymod_m=m)
   end subroutine new_molecular_grid_uniform_handymod

   !> Build a molecular grid with uniform sizes per atom
   !>
   !> The radial rule is midpoint quadrature on [-1, 1] followed by the
   !> HandyMod transform to [rmin, rmax]; the angular argument is an algebraic
   !> Lebedev degree
   !>
   !> @param[out] self            Initialised molecular grid
   !> @param[in]  mol             Molecular structure
   !> @param[in]  nrad            Number of midpoint radial nodes per atom
   !> @param[in]  lebedev_degree  Algebraic Lebedev angular degree
   !> @param[out] error           Propagated error
   !> @param[in]  rmin            Minimum HandyMod radius (bohr)
   !> @param[in]  rmax            Maximum HandyMod radius (bohr)
   !> @param[in]  m               Real-valued HandyMod transform parameter
   !> @param[in]  becke_k         Optional Becke polynomial stiffness (default 3)
   !> @param[in]  becke_ssf_a     Optional SSF cell parameter, replaces the polynomial
   !> @param[in]  becke_thr       Optional bare-Becke-weight prune threshold (default 0)
   !> @param[in]  arc_r           Optional per-shell angular rule: band-edge radii (bohr)
   !> @param[in]  arc_a           Optional per-shell angular rule: target arc spacing per band (bohr)
   !> @param[in]  nang_min        Optional floor on the per-shell Lebedev point count (default 110)
   !> @param[in]  nang_max        Optional cap on the per-shell Lebedev point count (default 5810)
   subroutine new_molecular_grid_uniform_qc_handymod(self, mol, nrad, lebedev_degree, error, &
         & rmin, rmax, m, becke_k, becke_ssf_a, becke_thr, arc_r, arc_a, nang_min, nang_max)
      !> Grid to initialise
      type(moist_math_grid_3d_molecular_type), intent(out) :: self
      !> Molecular structure
      type(structure_type), intent(in)  :: mol
      !> Number of radial midpoint nodes per atom
      integer, intent(in)  :: nrad
      !> Algebraic Lebedev angular degree
      integer, intent(in)  :: lebedev_degree
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Minimum HandyMod radius (bohr)
      real(wp), intent(in)  :: rmin
      !> Maximum HandyMod radius (bohr)
      real(wp), intent(in)  :: rmax
      !> Real-valued HandyMod transform parameter
      real(wp), intent(in)  :: m
      !> Optional Becke cutoff-polynomial stiffness `k` (default 3)
      !>
      !> Shared across the grid constructors, so absent must mean exactly
      !> Becke's k = 3
      integer, optional, intent(in)  :: becke_k
      !> Optional SSF cell-function parameter; replaces the Becke polynomial
      real(wp), optional, intent(in)  :: becke_ssf_a
      !> Optional bare-Becke-weight prune threshold
      !>
      !> Points whose partition weight is below this are dropped as owned
      !> by another atom; default 0 prunes nothing beyond the total-weight
      !> threshold `default_wthr`
      real(wp), optional, intent(in)  :: becke_thr
      !> Optional per-shell angular rule: ascending band-edge radii (bohr)
      real(wp), optional, intent(in)  :: arc_r(:)
      !> Optional per-shell angular rule: target arc spacing per band (bohr)
      !>
      !> One more entry than `arc_r`; when present it replaces the single
      !> per-atom order, and `lebedev_degree` is then unused for construction
      real(wp), optional, intent(in)  :: arc_a(:)
      !> Optional floor on the per-shell Lebedev point count (default 110)
      integer, optional, intent(in)  :: nang_min
      !> Optional cap on the per-shell Lebedev point count (default 5810)
      integer, optional, intent(in)  :: nang_max

      integer, allocatable :: nrad_atom(:), lebedev_degree_atom(:)

      if (nrad < 1) then
         call fatal_error(error, "molecular grid qc-HandyMod: nrad must be >= 1")
         return
      end if
      if (rmax <= rmin) then
         call fatal_error(error, "molecular grid qc-HandyMod: rmax must be greater than rmin")
         return
      end if
      if (m <= 0.0_wp) then
         call fatal_error(error, "molecular grid qc-HandyMod: m must be > 0")
         return
      end if

      allocate (nrad_atom(mol%nat), lebedev_degree_atom(mol%nat))
      nrad_atom = nrad
      lebedev_degree_atom = lebedev_degree

      call build_molecular_grid_qc_handymod(self, mol, nrad_atom, lebedev_degree_atom, &
         & default_wthr, error, rmin, rmax, m, becke_k=becke_k, becke_ssf_a=becke_ssf_a, &
         & becke_thr=becke_thr, arc_r=arc_r, arc_a=arc_a, nang_min=nang_min, nang_max=nang_max)
      if (allocated(error)) return
      call remember_molecular_settings(self, mol, recipe_qc_handymod, becke_k, becke_ssf_a, becke_thr, &
         & nrad=nrad, lebedev_degree=lebedev_degree, rmin=rmin, rmax=rmax, m=m, &
         & arc_r=arc_r, arc_a=arc_a, nang_min=nang_min, nang_max=nang_max)
   end subroutine new_molecular_grid_uniform_qc_handymod

   !> Internal worker: build atom-centered molecular grid
   !>
   !> Takes per-atom (nrad, nang) arrays and an optional radial clamp
   !>
   !> @param[out] self         Grid to initialise
   !> @param[in]  mol          Molecular structure
   !> @param[in]  nrad_atom    Per-atom radial sizes, shape (nat)
   !> @param[in]  nang_atom    Per-atom raw Lebedev point counts, shape (nat)
   !> @param[in]  wthr         Pruning threshold (bohr^3)
   !> @param[out] error        Propagated error
   !> @param[in]  rmin         Optional radial clamp (Chebyshev) or HandyMod lower bound (bohr)
   !> @param[in]  rmax         Optional radial clamp (Chebyshev) or HandyMod upper bound (bohr)
   !> @param[in]  handymod_m   Optional HandyMod exponent; absent keeps Chebyshev-2 behaviour
   !> @param[in]  becke_k      Optional Becke polynomial stiffness (default 3)
   !> @param[in]  becke_ssf_a  Optional SSF cell parameter, replaces the polynomial
   !> @param[in]  becke_thr    Optional bare-Becke-weight prune threshold (default 0)
   subroutine build_molecular_grid(self, mol, nrad_atom, nang_atom, wthr, &
         & error, rmin, rmax, handymod_m, becke_k, becke_ssf_a, becke_thr)
      !> Grid to initialise
      type(moist_math_grid_3d_molecular_type), intent(out) :: self
      !> Molecular structure
      type(structure_type), intent(in)  :: mol
      !> Per-atom radial sizes, shape (nat)
      integer, intent(in)  :: nrad_atom(:)
      !> Per-atom raw Lebedev point counts, shape (nat)
      integer, intent(in)  :: nang_atom(:)
      !> Pruning threshold (bohr^3)
      real(wp), intent(in)  :: wthr
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional radial clamps (Chebyshev) or HandyMod interval bounds (bohr)
      real(wp), optional, intent(in)  :: rmin, rmax
      !> Optional HandyMod exponent; absent keeps Chebyshev-2 behaviour
      integer, optional, intent(in)  :: handymod_m
      !> Optional Becke cutoff-polynomial stiffness (default 3)
      integer, optional, intent(in)  :: becke_k
      !> Optional SSF cell-function parameter; replaces the Becke polynomial
      real(wp), optional, intent(in)  :: becke_ssf_a
      !> Optional bare-Becke-weight prune threshold (default 0 = no pruning)
      real(wp), optional, intent(in)  :: becke_thr

      integer  :: nat, iat, iz, ir, il, ig, nr, nl, max_pts, order
      integer(int64) :: max_pts64
      integer, allocatable :: numbers(:)
      real(wp) :: p, r, full_weight, rlo, rhi
      integer :: hm
      logical :: use_handymod
      real(wp) :: bthr_val
      real(wp), allocatable :: radii(:), rad_w(:)
      real(wp), allocatable :: ang_xyz(:, :), ang_w(:)
      real(wp), allocatable :: tmp_xyz(:, :), tmp_w(:), tmp_bw(:)
      integer, allocatable :: tmp_atom(:)
      real(wp) :: bweights_buf(mol%nat)
      real(wp), parameter :: four_pi = 4.0_wp*pi

      nat = mol%nat
      use_handymod = present(handymod_m)
      rlo = 0.0_wp
      rhi = 0.0_wp
      hm = 0
      bthr_val = 0.0_wp
      if (present(becke_thr)) bthr_val = becke_thr

      if (use_handymod) then
         if (.not. present(rmin) .or. .not. present(rmax)) then
            call fatal_error(error, "molecular grid HandyMod: rmin and rmax are required")
            return
         end if
         rlo = rmin
         rhi = rmax
         hm = handymod_m
      end if

      ! Resolve atomic numbers once (mctc_io indexes through mol%id)
      allocate (numbers(nat))
      do iat = 1, nat
         numbers(iat) = mol%num(mol%id(iat))
      end do

      ! Upper bound on grid size before pruning, summed in int64 so a huge
      ! per-atom request overflows the check below rather than `max_pts`
      max_pts64 = 0_int64
      do iat = 1, nat
         max_pts64 = max_pts64 + int(nrad_atom(iat), int64)*int(nang_atom(iat), int64)
      end do
      call check_point_count(max_pts64, error)
      if (allocated(error)) return
      max_pts = int(max_pts64)

      allocate (tmp_xyz(3, max_pts), tmp_w(max_pts), tmp_bw(max_pts), tmp_atom(max_pts))

      ig = 0
      do iat = 1, nat
         iz = numbers(iat)
         nr = nrad_atom(iat)
         nl = nang_atom(iat)

         allocate (radii(nr), rad_w(nr))
         if (use_handymod) then
            call handymod_radii(nr, rlo, rhi, hm, radii, rad_w, error)
            if (allocated(error)) then
               deallocate (radii, rad_w)
               return
            end if
         else
            ! Chebyshev-2 radial scale p (Becke 1988 convention)
            if (iz == 1) then
               p = covalent_rad(1)
            else
               p = 0.5_wp*covalent_rad(iz)
            end if
            call chebyshev2_radii(nr, p, radii, rad_w)
         end if

         ! Lebedev angular grid (use existing umbrella getter)
         call lebedev_order_from_num(nl, order, error)
         if (allocated(error)) then
            deallocate (radii, rad_w)
            return
         end if
         allocate (ang_xyz(3, nl), ang_w(nl))
         call get_angular_grid(order, ang_xyz, ang_w, error)
         if (allocated(error)) then
            deallocate (radii, rad_w, ang_xyz, ang_w)
            return
         end if

         ! Combine radial x angular, shift to atom center, apply Becke weight
         do ir = 1, nr
            r = radii(ir)
            if (.not. use_handymod) then
               if (present(rmin)) then
                  if (r < rmin) cycle
               end if
               if (present(rmax)) then
                  if (r > rmax) cycle
               end if
            end if
            do il = 1, nl
               ig = ig + 1
               tmp_xyz(1, ig) = r*ang_xyz(1, il) + mol%xyz(1, iat)
               tmp_xyz(2, ig) = r*ang_xyz(2, il) + mol%xyz(2, iat)
               tmp_xyz(3, ig) = r*ang_xyz(3, il) + mol%xyz(3, iat)
               full_weight = rad_w(ir)*ang_w(il)*four_pi
               call becke_weights(tmp_xyz(:, ig), nat, mol%xyz, numbers, &
                  & bweights_buf, stiffness=becke_k, ssf_a=becke_ssf_a)
               tmp_bw(ig) = bweights_buf(iat)
               tmp_w(ig) = full_weight*bweights_buf(iat)
               tmp_atom(ig) = iat
            end do
         end do

         deallocate (radii, rad_w, ang_xyz, ang_w)
      end do

      call finalise_grid(self, nat, nrad_atom, nang_atom, &
         & ig, tmp_xyz, tmp_w, tmp_atom, wthr, raw_bw=tmp_bw, bthr=bthr_val)

      deallocate (tmp_xyz, tmp_w, tmp_bw, tmp_atom, numbers)
   end subroutine build_molecular_grid

   !> Internal worker: build a midpoint-HandyMod molecular grid
   !>
   !> @param[out] self                 Grid to initialise
   !> @param[in]  mol                  Molecular structure
   !> @param[in]  nrad_atom            Per-atom radial midpoint sizes, shape (nat)
   !> @param[in]  lebedev_degree_atom  Per-atom algebraic Lebedev degrees, shape (nat)
   !> @param[in]  wthr                 Pruning threshold (bohr^3)
   !> @param[out] error                Propagated error
   !> @param[in]  rmin                 Minimum HandyMod radius (bohr)
   !> @param[in]  rmax                 Maximum HandyMod radius (bohr)
   !> @param[in]  m                    Real-valued HandyMod transform parameter
   !> @param[in]  becke_k              Optional Becke polynomial stiffness (default 3)
   !> @param[in]  becke_ssf_a          Optional SSF cell parameter, replaces the polynomial
   !> @param[in]  becke_thr            Optional bare-Becke-weight prune threshold (default 0)
   !> @param[in]  arc_r                Optional per-shell angular rule: band-edge radii (bohr)
   !> @param[in]  arc_a                Optional per-shell angular rule: target arc spacing per band (bohr)
   !> @param[in]  nang_min             Optional floor on the per-shell Lebedev point count (default 110)
   !> @param[in]  nang_max             Optional cap on the per-shell Lebedev point count (default 5810)
   subroutine build_molecular_grid_qc_handymod(self, mol, nrad_atom, lebedev_degree_atom, &
         & wthr, error, rmin, rmax, m, becke_k, becke_ssf_a, becke_thr, &
         & arc_r, arc_a, nang_min, nang_max)
      !> Grid to initialise
      type(moist_math_grid_3d_molecular_type), intent(out) :: self
      !> Molecular structure
      type(structure_type), intent(in)  :: mol
      !> Per-atom radial midpoint sizes, shape (nat)
      integer, intent(in)  :: nrad_atom(:)
      !> Per-atom algebraic Lebedev degrees, shape (nat)
      integer, intent(in)  :: lebedev_degree_atom(:)
      !> Pruning threshold (bohr^3)
      real(wp), intent(in)  :: wthr
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Minimum HandyMod radius (bohr)
      real(wp), intent(in)  :: rmin
      !> Maximum HandyMod radius (bohr)
      real(wp), intent(in)  :: rmax
      !> Real-valued HandyMod transform parameter
      real(wp), intent(in)  :: m
      !> Optional Becke cutoff-polynomial stiffness (default 3)
      integer, optional, intent(in)  :: becke_k
      !> Optional SSF cell-function parameter; replaces the Becke polynomial
      real(wp), optional, intent(in)  :: becke_ssf_a
      !> Optional bare-Becke-weight prune threshold (default 0 = no pruning)
      real(wp), optional, intent(in)  :: becke_thr
      !> Optional per-shell angular rule: ascending band-edge radii (bohr)
      !>
      !> Present together with `arc_a` it replaces the single per-atom
      !> Lebedev order with one chosen shell by shell from the target arc
      !> spacing
      real(wp), optional, intent(in)  :: arc_r(:)
      !> Optional per-shell angular rule: target arc spacing per band (bohr)
      !>
      !> One more entry than `arc_r`
      real(wp), optional, intent(in)  :: arc_a(:)
      !> Optional floor on the per-shell Lebedev point count (default 110)
      integer, optional, intent(in)  :: nang_min
      !> Optional cap on the per-shell Lebedev point count (default 5810)
      integer, optional, intent(in)  :: nang_max

      integer :: nat, iat, ir, il, ig, nr, nl, max_pts, order, order_prev
      integer :: nlo, nhi
      integer(int64) :: max_pts64
      integer, allocatable :: numbers(:), nang_atom(:), shell_order(:)
      real(wp) :: r, full_weight, bthr_val
      logical :: use_pershell
      real(wp), allocatable :: radii(:), rad_w(:)
      real(wp), allocatable :: ang_xyz(:, :), ang_w(:)
      real(wp), allocatable :: tmp_xyz(:, :), tmp_w(:), tmp_bw(:)
      integer, allocatable :: tmp_atom(:)
      real(wp) :: bweights_buf(mol%nat)
      real(wp), parameter :: four_pi = 4.0_wp*pi

      nat = mol%nat
      bthr_val = 0.0_wp
      if (present(becke_thr)) bthr_val = becke_thr
      if (size(nrad_atom) < nat .or. size(lebedev_degree_atom) < nat) then
         call fatal_error(error, "molecular grid qc-HandyMod: per-atom size arrays are too small")
         return
      end if

      ! The per-shell rule is opt-in on arc_a; with it absent every line
      ! below reduces to the single-order-per-atom path unchanged
      use_pershell = present(arc_a)
      nlo = default_shell_nang_min
      nhi = default_shell_nang_max
      if (present(nang_min)) nlo = nang_min
      if (present(nang_max)) nhi = nang_max
      if (use_pershell) then
         if (.not. present(arc_r)) then
            call fatal_error(error, "molecular grid qc-HandyMod: arc_a requires arc_r")
            return
         end if
         call validate_arc_bands(arc_r, arc_a, nlo, nhi, error)
         if (allocated(error)) return
      end if

      allocate (numbers(nat), nang_atom(nat))
      do iat = 1, nat
         numbers(iat) = mol%num(mol%id(iat))
         call lebedev_degree_to_order(lebedev_degree_atom(iat), order, error)
         if (allocated(error)) return
         nang_atom(iat) = grid_size(order)
      end do

      ! Under the per-shell rule the point count is no longer
      ! nrad * nang_atom, so the raw-buffer bound has to be summed from the
      ! per-shell orders, costing one extra evaluation of the radial rule
      ! per atom and keeping the buffer exact rather than merely sufficient
      max_pts64 = 0_int64
      if (use_pershell) then
         do iat = 1, nat
            nr = nrad_atom(iat)
            allocate (radii(nr), rad_w(nr))
            call handymod_midpoint_radii(nr, rmin, rmax, m, radii, rad_w, error)
            if (allocated(error)) then
               deallocate (radii, rad_w)
               return
            end if
            do ir = 1, nr
               order = lebedev_order_from_arc(radii(ir), &
                  & arc_target_for_shell(radii(ir), arc_r, arc_a), nlo, nhi)
               max_pts64 = max_pts64 + int(grid_size(order), int64)
            end do
            deallocate (radii, rad_w)
         end do
      else
         do iat = 1, nat
            max_pts64 = max_pts64 + int(nrad_atom(iat), int64)*int(nang_atom(iat), int64)
         end do
      end if
      call check_point_count(max_pts64, error)
      if (allocated(error)) return
      max_pts = int(max_pts64)

      allocate (tmp_xyz(3, max_pts), tmp_w(max_pts), tmp_bw(max_pts), tmp_atom(max_pts))

      ig = 0
      do iat = 1, nat
         nr = nrad_atom(iat)
         call lebedev_degree_to_order(lebedev_degree_atom(iat), order, error)
         if (allocated(error)) return
         nl = grid_size(order)

         allocate (radii(nr), rad_w(nr), shell_order(nr))
         call handymod_midpoint_radii(nr, rmin, rmax, m, radii, rad_w, error)
         if (allocated(error)) then
            deallocate (radii, rad_w, shell_order)
            return
         end if

         if (use_pershell) then
            do ir = 1, nr
               shell_order(ir) = lebedev_order_from_arc(radii(ir), &
                  & arc_target_for_shell(radii(ir), arc_r, arc_a), nlo, nhi)
            end do
            ! Report the largest order actually used, so nang_per_atom stays
            ! the honest summary of the atom's angular resolution
            nang_atom(iat) = grid_size(maxval(shell_order(1:nr)))
         else
            shell_order(1:nr) = order
         end if

         ! The radial nodes ascend and the arc rule is monotone in r, so the
         ! order changes only a handful of times per atom; refetching the
         ! Lebedev grid on change alone keeps this to ~10 calls instead of one
         ! per shell
         order_prev = 0
         do ir = 1, nr
            r = radii(ir)
            if (shell_order(ir) /= order_prev) then
               if (allocated(ang_xyz)) deallocate (ang_xyz, ang_w)
               nl = grid_size(shell_order(ir))
               allocate (ang_xyz(3, nl), ang_w(nl))
               call get_angular_grid(shell_order(ir), ang_xyz, ang_w, error)
               if (allocated(error)) then
                  deallocate (radii, rad_w, shell_order, ang_xyz, ang_w)
                  return
               end if
               order_prev = shell_order(ir)
            end if
            do il = 1, nl
               ig = ig + 1
               tmp_xyz(1, ig) = mol%xyz(1, iat) + r*ang_xyz(1, il)
               tmp_xyz(2, ig) = mol%xyz(2, iat) + r*ang_xyz(2, il)
               tmp_xyz(3, ig) = mol%xyz(3, iat) + r*ang_xyz(3, il)
               full_weight = rad_w(ir)*ang_w(il)*four_pi
               call becke_weights(tmp_xyz(:, ig), nat, mol%xyz, numbers, bweights_buf, &
                  & stiffness=becke_k, ssf_a=becke_ssf_a)
               tmp_bw(ig) = bweights_buf(iat)
               tmp_w(ig) = full_weight*bweights_buf(iat)
               tmp_atom(ig) = iat
            end do
         end do

         deallocate (radii, rad_w, shell_order)
         if (allocated(ang_xyz)) deallocate (ang_xyz, ang_w)
      end do

      call finalise_grid(self, nat, nrad_atom, nang_atom, ig, tmp_xyz, tmp_w, tmp_atom, &
         & wthr, raw_bw=tmp_bw, bthr=bthr_val)

      deallocate (tmp_xyz, tmp_w, tmp_bw, tmp_atom, numbers, nang_atom)
   end subroutine build_molecular_grid_qc_handymod

   !> Map a raw Lebedev point count to its algebraic Lebedev degree
   !>
   !> The constructor selects its angular rule by algebraic degree, while
   !> the rest of moist (and every existing caller) speaks in raw point
   !> counts
   !> Inverts `grid_size` so a caller holding a point count can reach the
   !> degree-based constructor without duplicating the table
   !>
   !> @param[in]  num             Raw Lebedev point count (6, 14, ..., 5810)
   !> @param[out] lebedev_degree  Algebraic Lebedev degree of that grid
   !> @param[out] error           Propagated error for an unsupported count
   subroutine lebedev_degree_from_num(num, lebedev_degree, error)
      !> Raw Lebedev point count
      integer, intent(in) :: num
      !> Algebraic Lebedev degree
      integer, intent(out) :: lebedev_degree
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: order

      lebedev_degree = 0
      call lebedev_order_from_num(num, order, error)
      if (allocated(error)) return
      if (order < 1 .or. order > size(lebedev_degree_table)) then
         call fatal_error(error, "molecular grid: Lebedev count has no tabulated degree")
         return
      end if
      lebedev_degree = lebedev_degree_table(order)
   end subroutine lebedev_degree_from_num

   !> Map an algebraic Lebedev degree to the existing grid-order index
   !>
   !> @param[in]  lebedev_degree  Algebraic Lebedev angular degree
   !> @param[out] order           Existing moist Lebedev order index
   !> @param[out] error           Error handling
   subroutine lebedev_degree_to_order(lebedev_degree, order, error)
      !> Algebraic Lebedev angular degree
      integer, intent(in) :: lebedev_degree
      !> Existing moist Lebedev order index
      integer, intent(out) :: order
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i

      order = 0
      do i = 1, size(lebedev_degree_table)
         if (lebedev_degree_table(i) == lebedev_degree) then
            order = i
            return
         end if
      end do

      call fatal_error(error, "molecular grid qc-HandyMod: unsupported Lebedev degree")
   end subroutine lebedev_degree_to_order

   !> Guard a raw-buffer point-count accumulator against integer overflow
   !>
   !> `total` is computed in int64 by the caller, which sums per-atom/per-shell
   !> point counts before this is called, so the overflow is caught before the
   !> subsequent `allocate(tmp_xyz(3, max_pts), ...)` at default-integer size
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

   !> Compact the raw grid buffer into the final `moist_math_grid_3d_molecular_type`
   !>
   !> Drops points with |w| < wthr and rebuilds `atom_offset`
   !>
   !> @param[out] self       Grid to fill
   !> @param[in]  nat        Number of atoms
   !> @param[in]  nrad_atom  Per-atom radial sizes, shape (nat)
   !> @param[in]  nang_atom  Per-atom raw Lebedev point counts, shape (nat)
   !> @param[in]  nraw       Number of raw (pre-pruning) points
   !> @param[in]  raw_xyz    Raw point coordinates, shape (3, nraw), bohr
   !> @param[in]  raw_w      Raw quadrature weights, shape (nraw), bohr^3
   !> @param[in]  raw_atom   Raw point owning atom, shape (nraw), one based
   !> @param[in]  wthr       Pruning threshold on |raw_w| (bohr^3)
   !> @param[in]  raw_bw     Optional bare Becke partition weight of each raw point, shape (nraw)
   !> @param[in]  bthr       Optional Becke-weight prune threshold
   pure subroutine finalise_grid(self, nat, nrad_atom, nang_atom, &
         & nraw, raw_xyz, raw_w, raw_atom, wthr, raw_bw, bthr)
      !> Grid to fill
      type(moist_math_grid_3d_molecular_type), intent(out) :: self
      !> Number of atoms
      integer, intent(in) :: nat
      !> Per-atom radial sizes, shape (nat)
      integer, intent(in) :: nrad_atom(:)
      !> Per-atom raw Lebedev point counts, shape (nat)
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
      !> Bare Becke partition weight of each raw point, shape (nraw)
      !>
      !> Dimensionless `w_A(r)` BEFORE the `r^2 dr dOmega` Jacobian, which
      !> is what makes it usable as a redundancy threshold: sums to 1 over
      !> atoms and does not vanish near a nucleus the way the total weight
      !> does
      real(wp), intent(in), optional :: raw_bw(:)
      !> Becke-weight prune threshold
      !>
      !> Points another atom owns more strongly than this are dropped; zero
      !> or absent prunes nothing beyond `wthr`
      real(wp), intent(in), optional :: bthr

      integer :: i, iat, ngrid
      logical :: prune_bw
      logical, allocatable :: keep(:)

      ! Written as an explicit opt-in rather than as raw_bw >= bthr so that
      ! the default path is bit-identical to the unpruned build for every
      ! possible floating-point value of the Becke weight
      prune_bw = .false.
      if (present(raw_bw) .and. present(bthr)) then
         if (bthr > 0.0_wp) prune_bw = .true.
      end if

      allocate (keep(nraw))
      if (prune_bw) then
         do i = 1, nraw
            keep(i) = abs(raw_w(i)) >= wthr .and. raw_bw(i) >= bthr
         end do
      else
         do i = 1, nraw
            keep(i) = abs(raw_w(i)) >= wthr
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
      ! raw_atom is monotonically non-decreasing because we emit atom-by-atom
      ! Default every entry to the end sentinel so atoms with zero points
      ! get a well-defined (empty) range; overwritten as points arrive
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

   !> Free all component arrays, idempotent
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
      self%ngrid = 0
      self%npts_k = 0
      self%nkx = 0; self%nky = 0; self%nkz = 0
      self%dkx = 0.0_wp; self%dky = 0.0_wp; self%dkz = 0.0_wp
      self%kref = 0.0_wp
      self%has_kgrid = .false.
      self%centered_period = .false.
      self%kref_offset = 0.0_wp
      ! Invalidate any trafo still bound to this grid from before the destroy
      self%geom_generation = self%geom_generation + 1
   end subroutine molecular_grid_destroy

   !> Configure the uniform reciprocal (k) grid used by the NUFFT
   !>
   !> - Single-resolution interface: the caller supplies only the real-space
   !>   spacing `dr` of the implied uniform box (this fixes the Nyquist
   !>   reach `k_max = pi/dr`), and the box is auto-sized to enclose every
   !>   grid point
   !> - For each axis the box length is `L = extent + 2*buffer`, the mode
   !>   count `N` is the smallest even integer with `N*dr >= L`, and the
   !>   spacing is `dk = 2*pi/(N*dr)`
   !> - Modes follow FINUFFT's CMCL layout (kx fastest, frequencies
   !>   -N/2 .. N/2-1)
   !> - Without an explicit `reference` the phase reference defaults to
   !>   `point(1)`, matching the Cartesian grid and `potential_3d`
   !> - `realize_molecular` passes the solute centroid as `reference` for
   !>   the auto-configured domain path (`auto_kgrid`, the type default); a
   !>   grid whose reciprocal grid was set explicitly (direct constructors)
   !>   keeps whichever convention that call chose across every later
   !>   update/rebuild
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
      !> Optional requested FINUFFT relative tolerance; largest cost lever on this backend
      real(wp), optional, intent(in)    :: nufft_tol
      !> Optional phase reference; sizes a symmetric period about this point
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
         ! Check before CEILING: an overlarge count can overflow its integer
         ! result and otherwise become the two-mode minimum below
         if (box/real(max_kgrid_modes_per_axis, wp) > dr) then
            call fatal_error(error, "molecular grid: auto-sized reciprocal grid is too "// &
               & "large (point cloud spans too wide for this dr); clamp the radial "// &
               & "extent (rmax) when building the grid or use a coarser dr")
            return
         end if
         ! Smallest even mode count whose implied box N*dr covers the extent
         nk(d) = ceiling(box/dr)
         nk(d) = max(nk(d), 2)
         if (mod(nk(d), 2) /= 0) nk(d) = nk(d) + 1
      end do

      if (any(nk > max_kgrid_modes_per_axis)) then
         call fatal_error(error, "molecular grid: auto-sized reciprocal grid is too "// &
            & "large (point cloud spans too wide for this dr); clamp the radial "// &
            & "extent (rmax) when building the grid or use a coarser dr")
         return
      end if

      self%nkx = nk(1); self%nky = nk(2); self%nkz = nk(3)
      self%dkx = two_pi/(real(nk(1), wp)*dr)
      self%dky = two_pi/(real(nk(2), wp)*dr)
      self%dkz = two_pi/(real(nk(3), wp)*dr)
      self%npts_k = nk(1)*nk(2)*nk(3)
      ! Default phase reference = point(1) (matches the Cartesian/potential_3d
      ! DFT) unless an explicit reference overrides it below
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
      ! Reconfiguring the reciprocal grid invalidates any trafo already
      ! prepared against the previous nkx/nky/nkz, dk* or kref
      self%geom_generation = self%geom_generation + 1
   end subroutine molecular_grid_set_kgrid

   !> Reciprocal-space coordinate (1/bohr) of k-point j (1..npts_k)
   !>
   !> Inverts the CMCL flat layout `j-1 = ax + nkx*(ay + nky*az)` (kx
   !> fastest) and maps each axis index to its signed frequency
   !> `m = a - N/2` (modes run -N/2 .. N/2-1), returning
   !> `(m_x*dkx, m_y*dky, m_z*dkz)`
   !> Requires the reciprocal grid to have been configured by
   !> `molecular_grid_set_kgrid`
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

   !> Allocate a `moist_math_grid_3d_molecular_trafo_type` bound to `self`
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

   !> Initialize a transform engine bound to `grid` (binds the pointer only)
   !>
   !> Call `trafo%prepare(nv)` afterwards to build the FINUFFT plans, once
   !> the grid's reciprocal grid has been set with `molecular_grid_set_kgrid`
   !> The grid must outlive the trafo
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

   !> Build the FINUFFT guru plans for `ntrans = nv` simultaneous transforms
   !>
   !> Two mathematically identical routes are available; both use the shifted
   !> real-space coordinates `x_j = xyz(:,j) - kref` so the phase reference
   !> matches the Cartesian FFT and `potential_3d`
   !>
   !> * **Type 1/2 (default, `want_type12`).**  The k-points are not an
   !>   arbitrary cloud: `molecular_grid_set_kgrid` builds them as the FFT dual
   !>   `k = m*dk` of a uniform box of side `nk*dr`, with `dk = 2*pi/(nk*dr)`
   !>   per axis and integer modes `m = -N/2 .. N/2-1`.  Writing the phase as
   !>   `k.x = m*(dk*x)` turns the transform into FINUFFT's type-1 (nonuniform
   !>   points -> uniform modes, forward) and type-2 (uniform modes ->
   !>   nonuniform points, backward) pair on the scaled coordinates
   !>   `dk*x`.  This is an identity, not an approximation: `exp(i*m*x)` with
   !>   integer `m` is exactly `2*pi`-periodic in `x`, so the periodic type-1/2
   !>   definition reproduces the non-periodic type-3 sum for any point spread
   !>   Type 1/2 sizes its internal fine grid from the mode count alone, while
   !>   type 3 sizes it from the spatial-half-width x frequency-half-width
   !>   product and additionally prephases and deconvolves, so the lattice
   !>   route is the cheaper of the two
   !>   FINUFFT's CMCL mode layout (`modeord = 0`, x fastest, `-N/2 ..
   !>   N/2-1`) is exactly the layout `molecular_grid_kpoint` and the
   !>   k-point loop below assume, so the flat ordering of `f_k` is unchanged
   !> * **Type 3 (`want_type12 = .false.`).**  The reference route: explicit
   !>   nonuniform k coordinates, forward `iflag = -1`, backward `iflag = +1`
   !>   Kept so the fast path can be A/B-compared elementwise in one process,
   !>   and used automatically if the scaled coordinates would leave FINUFFT's
   !>   portable `[-3pi, 3pi]` input window (they cannot, given
   !>   `nk*dr >= extent`, but the check is free and the fallback is safe)
   !>
   !> The grid's reciprocal grid must already be configured by
   !> `molecular_grid_set_kgrid`.  Re-preparing tears down any existing plans
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

      ! Idempotent: drop any plans/scratch from a previous prepare
      call molecular_trafo_free_plans(self)

      ngrid = self%grid%ngrid

      ! Shifted real-space coordinates
      ! The shift to kref (point(1) unless molecular_grid_set_kgrid was
      ! given an explicit reference) keeps the convention identical to the
      ! Cartesian FFT on either route
      npts_k = self%grid%npts_k
      allocate (self%xj(ngrid), self%yj(ngrid), self%zj(ngrid))
      do j = 1, ngrid
         self%xj(j) = self%grid%xyz(1, j) - self%grid%kref(1)
         self%yj(j) = self%grid%xyz(2, j) - self%grid%kref(2)
         self%zj(j) = self%grid%xyz(3, j) - self%grid%kref(3)
      end do

      ! Decide the transform family
      ! On the type-1/2 route the coordinates are scaled by the per-axis dk
      ! so k.x becomes m.(dk*x); the scaled spread is 2*pi*extent/(nk*dr)
      ! <= 2*pi by construction, so the fallback below never fires for a
      ! grid built by set_kgrid
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
      ! CMCL mode ordering (-N/2 .. N/2-1, x fastest): the layout the type-1/2
      ! route depends on, and the one molecular_grid_kpoint inverts
      opts%modeord = 0
      ! Explicit worker count, the Cartesian policy; left at 0 FINUFFT would
      ! take OMP_NUM_THREADS or the core count and ignore omp_set_num_threads
      opts%nthreads = 1_c_int
!$    if (.not. omp_in_parallel()) opts%nthreads = int(omp_get_max_threads(), c_int)

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

      ! Only now, with both plans fully built, does the trafo become usable;
      ! `ntrans` and `geom_generation` gate `check_ready`, so setting them
      ! here (rather than up front) keeps a partially-built engine unusable
      self%ntrans = ntrans
      self%geom_generation = self%grid%geom_generation
   end subroutine molecular_trafo_prepare

   !> Tear down the FINUFFT plans and free coordinate/strength scratch
   !>
   !> Idempotent module helper shared by `prepare` (re-prepare) and `destroy`
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

   !> Release the plans/scratch and detach the trafo from its grid, idempotent
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
   !> Real block `f_r(ngrid, nv)` -> reciprocal block `f_k(npts_k, nv)`, with
   !> strengths `w_j*f(r_j)` so the quadrature weights are baked in just like
   !> the Cartesian `dV*FFT`.  Requires `prepare(nv)` first; `nv` must equal
   !> the prepared `ntrans`, and a violated precondition or a FINUFFT
   !> failure is reported through `error`
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
   !> Reciprocal block `f_k(npts_k, nv)` -> real block `f_r(ngrid, nv)`,
   !> scaled by `dkx*dky*dkz/(2*pi)^3` to match the continuous inverse-FT
   !> normalisation (the Cartesian `1/Vbox`)
   !> Requires `prepare(nv)` first; a violated precondition or a FINUFFT
   !> failure is reported through `error`
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
      ! Argument order follows FINUFFT's execute(plan, cj, fk): on the type-2
      ! route cj is the (nonuniform) OUTPUT and f_k the mode INPUT, while
      ! on the type-3 route the k-space block is the source ("cj" slot) and
      ! the real-space values the target ("fk" slot); both write real-space
      ! values into self%cj
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
   !> The trafo must be prepared, bound to the same grid geometry it was
   !> prepared against, and the caller's block sizes must match that grid;
   !> reports a violated precondition through `error`
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

   !> Smallest admissible Lebedev order whose arc spacing at radius r meets a target spacing
   !>
   !> The arc spacing of an `n`-point Lebedev grid on a sphere of radius `r` is
   !> `r*sqrt(4*pi/n)`, so meeting `arc_target` needs
   !> `n >= 4*pi*r**2/arc_target**2`.  Sizes with negative quadrature weights
   !> are skipped and the result is clamped to `[nang_min, nang_max]`; when the
   !> target cannot be met within the cap the largest admissible size is
   !> returned, so the caller must read the ACHIEVED spacing back rather than
   !> assume the target was honoured
   !>
   !> @param[in] r           Shell radius (bohr)
   !> @param[in] arc_target  Requested arc spacing (bohr), must be positive
   !> @param[in] nang_min    Floor on the returned point count
   !> @param[in] nang_max    Cap on the returned point count
   pure function lebedev_order_from_arc(r, arc_target, nang_min, nang_max) result(order)
      !> Shell radius in bohr
      real(wp), intent(in) :: r
      !> Requested arc spacing in bohr
      real(wp), intent(in) :: arc_target
      !> Floor on the point count
      integer, intent(in) :: nang_min
      !> Cap on the point count
      integer, intent(in) :: nang_max
      !> Order index into `grid_size`, 0 if no admissible size exists
      integer :: order

      !> Loop index over the supported Lebedev sizes
      integer :: i
      !> Smallest point count that meets the target at this radius
      real(wp) :: needed

      order = 0
      if (arc_target <= 0.0_wp) return
      needed = 4.0_wp*pi*r*r/(arc_target*arc_target)

      ! `grid_size` is ascending, so the first admissible entry at or above the
      ! requirement is the smallest one
      do i = 1, size(grid_size)
         if (any(grid_size(i) == lebedev_negative_weight_sizes)) cycle
         if (grid_size(i) < nang_min) cycle
         if (grid_size(i) > nang_max) exit
         if (real(grid_size(i), wp) >= needed) then
            order = i
            return
         end if
      end do

      ! Target unreachable inside the cap: fall back to the largest admissible
      ! size and let the caller observe the shortfall
      do i = size(grid_size), 1, -1
         if (any(grid_size(i) == lebedev_negative_weight_sizes)) cycle
         if (grid_size(i) > nang_max) cycle
         if (grid_size(i) < nang_min) cycle
         order = i
         return
      end do
   end function lebedev_order_from_arc

   !> Target arc spacing for a shell, from a piecewise-constant band rule
   !>
   !> `arc_r` holds ascending band-edge radii and `arc_a` the target spacing
   !> in each band, with one more entry than `arc_r` (the last applies
   !> beyond the final edge)
   !> A shell at exactly an edge belongs to the inner band
   !>
   !> @param[in] r      Shell radius (bohr)
   !> @param[in] arc_r  Ascending band-edge radii (bohr)
   !> @param[in] arc_a  Target arc spacing per band (bohr), size(arc_r) + 1
   pure function arc_target_for_shell(r, arc_r, arc_a) result(a)
      !> Shell radius in bohr
      real(wp), intent(in) :: r
      !> Ascending band-edge radii in bohr
      real(wp), intent(in) :: arc_r(:)
      !> Target arc spacing per band in bohr
      real(wp), intent(in) :: arc_a(:)
      !> Target arc spacing at `r`
      real(wp) :: a

      !> Loop index over the band edges
      integer :: i

      a = arc_a(size(arc_a))
      do i = 1, size(arc_r)
         if (r <= arc_r(i)) then
            a = arc_a(i)
            return
         end if
      end do
   end function arc_target_for_shell

   !> Validate a per-shell angular band rule
   !>
   !> @param[in]  arc_r  Band-edge radii (bohr); must be ascending and positive
   !> @param[in]  arc_a  Target spacings (bohr); must be positive
   !> @param[in]  nlo    Floor on the point count
   !> @param[in]  nhi    Cap on the point count
   !> @param[out] error  Set when the rule is not usable
   subroutine validate_arc_bands(arc_r, arc_a, nlo, nhi, error)
      !> Band-edge radii in bohr
      real(wp), intent(in) :: arc_r(:)
      !> Target spacings in bohr
      real(wp), intent(in) :: arc_a(:)
      !> Floor on the point count
      integer, intent(in) :: nlo
      !> Cap on the point count
      integer, intent(in) :: nhi
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Loop index over the band edges
      integer :: i

      if (size(arc_a) /= size(arc_r) + 1) then
         call fatal_error(error, "molecular grid per-shell angular rule: "// &
            & "arc_a must have exactly one more entry than arc_r")
         return
      end if
      if (any(arc_a <= 0.0_wp)) then
         call fatal_error(error, "molecular grid per-shell angular rule: "// &
            & "arc_a entries must be positive")
         return
      end if
      do i = 1, size(arc_r)
         if (arc_r(i) <= 0.0_wp) then
            call fatal_error(error, "molecular grid per-shell angular rule: "// &
               & "arc_r entries must be positive")
            return
         end if
         if (i > 1) then
            if (arc_r(i) <= arc_r(i - 1)) then
               call fatal_error(error, "molecular grid per-shell angular rule: "// &
                  & "arc_r entries must be strictly ascending")
               return
            end if
         end if
      end do
      if (nlo > nhi) then
         call fatal_error(error, "molecular grid per-shell angular rule: "// &
            & "nang_min exceeds nang_max")
         return
      end if
      if (lebedev_order_from_arc(1.0_wp, 1.0_wp, nlo, nhi) == 0) then
         call fatal_error(error, "molecular grid per-shell angular rule: "// &
            & "no supported Lebedev size lies within [nang_min, nang_max]")
         return
      end if
   end subroutine validate_arc_bands

end module moist_math_grid_3d_molecular
