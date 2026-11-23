module moist_cavity_drop_parameters
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use moist_model_parameters, only: moist_model_parameters_type
   use, intrinsic :: iso_fortran_env, only: output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_cavity_bornzeta, only: get_born_zeta
   use moist_cavity_drop_projection_batch, only: default_batch_size
   use moist_cavity_drop_projection_parameters, only: drop_seed_parameters_type, drop_solve_parameters_type, &
                                                     drop_polish_parameters_type, drop_finalize_parameters_type, &
                                                     drop_branch_parameters_type

   implicit none(type, external)

   public :: moist_cavity_drop_parameters_type
   public :: is_valid_proj_level
   public :: drop_seed_parameters_type, drop_branch_parameters_type
   public :: drop_solve_parameters_type, drop_polish_parameters_type, drop_finalize_parameters_type

   private

   !> Length of a `proj_level_label` result
   integer, parameter :: proj_level_label_len = 32

   !> Parameter container for DROP
   !> Extends the base parameter type with DROP-specific settings for
   !> grid generation, blending, barrier potentials, and optimization
   type, extends(moist_model_parameters_type) :: moist_cavity_drop_parameters_type

      !* ----------------------------- Grid Discretization ---------------------------- *!

      !> Number of Lebedev quadrature points per atomic sphere
      integer :: num_leb = 194

      !* ----------------------------- Optional properties ---------------------------- *!

      !> Compute all optional cavity properties below
      logical :: do_fine = .false.

      !> Compute the local grid point density
      logical :: do_grid_density = .false.
      !> Compute principal curvatures (and mean and Gaussian)
      logical :: do_curvature = .false.
      !> Store surface normal vectors at grid points
      logical :: do_normal = .false.
      !> Atomic areas via spatial partition of atomic contributions of the LSF
      logical :: do_atomarea = .false.

      !* ---------------------------------- Tolerance --------------------------------- *!

      !> Main tolerance controlling all numerical thresholds
      !>
      !> Derived tolerances (tightest to loosest):
      !>  - finalize%wleb_cut    = tolerance * 0.05 (quadrature weight cutoff)
      !>  - screening_threshold = tolerance * 0.1  (LSF screening; passed into lsf_model by DROP constructor)
      !>  - solve%tol           = tolerance        (projection convergence)
      !>  - branch%sep_cut      = tolerance * 10   (branch degeneracy)
      real(wp) :: tolerance = 1.0E-10_wp

      !> LSF screening threshold, derived from `tolerance` in
      !> `compute_drop_derived`
      !>
      !> - the DROP constructor pushes this value into the polymorphic
      !>   `lsf_model`, so the LSF's own screening caches stay consistent with
      !>   the cavity tolerance
      real(wp) :: screening_threshold = 1.0E-11_wp

      !* ------------------------- Objective Function Weights ------------------------- *!

      !> Weight w_a for anchor term
      real(wp) :: phi_alpha = 0.5_wp

      !* --------------------------- Projection Optimization -------------------------- *!

      !> Projection level:
      !> 1 = local batched KKT Newton only
      !> 3 = every anchor, coarse octree seeds (`seed%init_coarse`)
      !> 4 = every anchor, fine octree seeds (`seed%init_fine`)
      !> Any other value is refused
      !> A failed local solve is retried from the level's octree, the coarse one at level 1
      !> The coarse octree drops a leaf without sign change; the fine one probes its children once
      integer :: proj_level = 3

      !* ----------------------------- Batched projection ----------------------------- *!

      !> Anchors per projection chunk and flush target of the seeded pass (points)
      integer :: batch_size = default_batch_size

      !* ------------------------------ Projection stages ----------------------------- *!

      !> Seed step (octree surface coverage); leaf size and refinement derived from `proj_level`
      type(drop_seed_parameters_type) :: seed
      !> Solve step (batched KKT Newton); `tol` derived from `tolerance`
      type(drop_solve_parameters_type) :: solve
      !> Optimality polish
      type(drop_polish_parameters_type) :: polish
      !> Per-anchor finalize; `wleb_cut` derived from `tolerance`,
      !> the prune switch from `wleb_prune_level`
      type(drop_finalize_parameters_type) :: finalize

      !* ---------------------------------- Screening --------------------------------- *!

      !> Distance cutoff for grid point adj. list, derived as `4*rho_grid_h`
      real(wp) :: adj_list_grid_cutoff = 1.0_wp
      !> Below this atom count, the level set's candidate cell list collapses to a
      !> single full-scan cell
      !>
      !> - for small systems the per-cell fan-out reduces to "every atom in every
      !>   cell" anyway, so the build work is skipped and every query returns the
      !>   full atom list directly
      !> - handed to the level set through `set_full_scan_below`
      integer :: cell_grid_full_scan_below = 200

      !* ------------------------ Hard-sphere reference cavity ------------------------ *!

      !> Active iSwiG Gaussian zeta fitted to the analytical Born energy
      !> for the selected Lebedev grid size
      real(wp) :: iswig_xi_born = 0.0_wp

      !* ----------------------------- Switching Functions ---------------------------- *!

      !> Off-edge of critical level set weight switching transition
      real(wp) :: w_0ls_off = 0.1_wp
      !> On-edge of critical level set weight switching transition
      real(wp) :: w_0ls_on = 0.6_wp

      !> Off-edge of focal/branching weight switching transition, in the focal
      !> margin `beta_min/phi_alpha` (one on a flat surface)
      real(wp) :: w_0tra_off = 0.2_wp
      !> On-edge of focal/branching weight switching transition
      real(wp) :: w_0tra_on = 0.6_wp

      !> Off-edge of horizon weight switching transition, in the signed `n_surf . n_sph`
      !>
      !> - zero: closed exactly where the area scaling, proportional to the
      !>   cosine, vanishes
      real(wp) :: w_0hor_off = 0.0_wp
      !> On-edge of horizon weight switching transition
      !>
      !> - zero switches the horizon factor off
      real(wp) :: w_0hor_on = 0.15_wp
      !> Horizon band in use, derived: the explicit one, else a band below
      !> every cosine when switched off
      real(wp) :: hor_off = 0.0_wp
      real(wp) :: hor_on = 0.15_wp

      !> Function parameters
      real(wp) :: w_0ls_p = 0.8_wp
      real(wp) :: w_0ls_a = 1.6_wp

      !* ----------------------------- Grid point density ----------------------------- *!

      !> Grid density kernel length
      real(wp) :: rho_grid_h = 1.0_wp

      !* ---------------------------------- Branching --------------------------------- *!

      !> Branch window and reach; the window also bounds the branch search
      type(drop_branch_parameters_type) :: branch

      !* ------------------------------ Weight switching ------------------------------ *!

      !> Smooth switching on final Lebedev weights before the branch filter; the
      !> switch edges in `finalize` are derived from the level
      !>
      !> - Level 0: disabled (default)
      !> - Level 1: off at 1E-12, on at 1E-10
      !> - Level 2: off at 1E-10, on at 1E-8
      !> - Level 3: off at 1E-8,  on at 1E-6
      !> - Level 4: off at 1E-6,  on at 1E-4
      !> - Level 5: off at 1E-4,  on at 1E-2
      !> - Level 6: off at 1E-2,  on at 1E0
      integer :: wleb_prune_level = 0

   contains
      !> Initialize parameters to compiled defaults
      procedure :: init_defaults => init_cavity_drop_defaults
      !> Initialize parameters from constructor inputs
      procedure :: new => new_moist_cavity_drop_parameters_type
      !> Register parameters for JSON/TOML configuration parsing
      procedure :: register_entries => register_cavity_drop_entries
      !> Recompute all derived parameters from user-facing fields
      procedure :: compute_derived => compute_drop_derived
      !> Validate input and recompute derived values after file input
      procedure :: validate => compute_drop_derived
      !> Load parameters from a JSON or TOML file and recompute derived values
      procedure :: load_file => load_drop_file
      !> Select the fitted Born zeta value for the active Lebedev grid
      procedure :: select_born_zeta => select_born_zeta
      !> Print current parameter values to output
      procedure :: print => print_parameters
      !> Return a human-readable description for the current projection level
      procedure :: proj_level_label => get_proj_level_label
   end type moist_cavity_drop_parameters_type

contains

   !> Construct a new DROP parameters instance
   !>
   !> @param[inout] self Initialized parameter container
   !> @param[in]    nleb Number of Lebedev quadrature points
   !> @param[in]    tolerance main tolerance (derives solve tol, wleb_cut, branch_sep)
   !> @param[in]    proj_maxiter Maximum projection iterations
   !> @param[in]    proj_level Projection refinement level
   !> @param[in]    branch Branch window and reach
   !> @param[in]    rho_grid_h Grid density kernel length
   !> @param[in]    wleb_prune_level Weight switching level (0=off, 1-6=increasing aggressiveness)
   !> @param[in]    w_0hor_off Horizon switch off-edge
   !> @param[in]    w_0hor_on Horizon switch on-edge; zero switches it off
   !> @param[out]   error Invalid setting, or no fitted Born zeta for `nleb`
   subroutine new_moist_cavity_drop_parameters_type(self, &
                                                    nleb, tolerance, proj_maxiter, proj_level, &
                                                    branch, rho_grid_h, wleb_prune_level, &
                                                    w_0hor_off, w_0hor_on, &
                                                    error)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      integer, intent(in), optional :: nleb
      real(wp), intent(in), optional :: tolerance
      integer, intent(in), optional :: proj_maxiter
      integer, intent(in), optional :: proj_level
      type(drop_branch_parameters_type), intent(in), optional :: branch
      real(wp), intent(in), optional :: rho_grid_h
      integer, intent(in), optional :: wleb_prune_level
      real(wp), intent(in), optional :: w_0hor_off, w_0hor_on
      type(error_type), allocatable, intent(out) :: error

      ! Phase 1: Reset all fields to compiled defaults
      call self%init_defaults()

      ! Phase 2: Apply caller overrides
      if (present(nleb)) self%num_leb = nleb

      if (present(tolerance)) self%tolerance = tolerance
      if (present(proj_maxiter)) self%solve%max_iter = proj_maxiter
      if (present(proj_level)) self%proj_level = proj_level

      if (present(branch)) self%branch = branch
      if (present(rho_grid_h)) self%rho_grid_h = rho_grid_h
      if (present(wleb_prune_level)) self%wleb_prune_level = wleb_prune_level
      if (present(w_0hor_off)) self%w_0hor_off = w_0hor_off
      if (present(w_0hor_on)) self%w_0hor_on = w_0hor_on

      ! Phase 3: Compute derived parameters from final field values
      call self%compute_derived(error)
      if (allocated(error)) return

   end subroutine new_moist_cavity_drop_parameters_type

   !> Reset all user-facing parameters to compiled default values
   !>
   !> @param[inout] self Parameter container to reset
   subroutine init_cavity_drop_defaults(self)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      !> Fresh instance carrying the compiled defaults from the type declaration
      type(moist_cavity_drop_parameters_type) :: fresh

      ! Grid
      self%num_leb = fresh%num_leb
      self%do_fine = fresh%do_fine
      ! Optional properties
      self%do_grid_density = fresh%do_grid_density
      self%do_curvature = fresh%do_curvature
      self%do_normal = fresh%do_normal
      self%do_atomarea = fresh%do_atomarea
      ! Tolerance (master; solve%tol, finalize%wleb_cut, branch%sep_cut
      ! are derived from this in compute_derived; LSF screening lives on the LSF concrete)
      self%tolerance = fresh%tolerance
      ! Objective
      self%phi_alpha = fresh%phi_alpha
      ! Projection
      self%proj_level = fresh%proj_level
      ! Batched projection
      self%batch_size = fresh%batch_size
      ! Projection stages, derived tolerances included
      self%seed = fresh%seed
      self%solve = fresh%solve
      self%polish = fresh%polish
      self%finalize = fresh%finalize
      ! Screening
      self%cell_grid_full_scan_below = fresh%cell_grid_full_scan_below
      ! Switching
      self%w_0ls_off = fresh%w_0ls_off
      self%w_0ls_on = fresh%w_0ls_on
      self%w_0ls_p = fresh%w_0ls_p
      self%w_0ls_a = fresh%w_0ls_a
      self%w_0tra_off = fresh%w_0tra_off
      self%w_0tra_on = fresh%w_0tra_on
      self%w_0hor_off = fresh%w_0hor_off
      self%w_0hor_on = fresh%w_0hor_on
      self%hor_off = fresh%hor_off
      self%hor_on = fresh%hor_on
      self%wleb_prune_level = fresh%wleb_prune_level
      ! Density
      self%rho_grid_h = fresh%rho_grid_h
      ! Branching
      self%branch = fresh%branch
      !> Reset derived values; construction and file input recompute them
      self%screening_threshold = fresh%screening_threshold
      self%adj_list_grid_cutoff = fresh%adj_list_grid_cutoff
      self%iswig_xi_born = fresh%iswig_xi_born
   end subroutine init_cavity_drop_defaults

   !> Compute all derived parameters from the current user-facing fields
   !>
   !> Must be called after any modification of user-facing parameters
   !> (constructor or file load) to keep derived values consistent
   !>
   !> @param[inout] self Parameter container
   !> @param[out]   error Error if no fitted Born zeta exists for `num_leb`
   subroutine compute_drop_derived(self, error)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      !> Value under test is finite and within its bound
      logical :: valid

      valid = ieee_is_finite(self%tolerance)
      if (valid) valid = self%tolerance > 0.0_wp
      if (.not. valid) then
         call fatal_error(error, "Invalid DROP numerical parameters: tolerance must be finite and positive")
         return
      end if
      valid = ieee_is_finite(self%rho_grid_h)
      if (valid) valid = self%rho_grid_h > 0.0_wp
      if (.not. valid) then
         call fatal_error(error, "Invalid DROP numerical parameters: rho_grid_h must be finite and positive")
         return
      end if
      if (self%solve%max_iter < 1) then
         call fatal_error(error, "Invalid DROP numerical parameters: projection maxiter must be positive")
         return
      end if
      if (.not. is_valid_proj_level(self%proj_level)) then
         call fatal_error(error, "Invalid DROP numerical parameters: projection level must be 1, 3 or 4")
         return
      end if
      !> Octree leaf size and refinement of the level; the coarse octree is also
      !> the fallback of level 1
      if (self%proj_level == 4) then
         call self%seed%init_fine()
      else
         call self%seed%init_coarse()
      end if

      if (self%seed%max_boxes < 1) then
         call fatal_error(error, "Octree box budget must be positive")
         return
      end if
      if (self%seed%max_survivors < 1) then
         call fatal_error(error, "Octree survivor budget must be positive")
         return
      end if
      if (self%seed%max_depth < 1 .or. self%seed%max_depth > self%seed%max_addressable_depth) then
         call fatal_error(error, "Octree depth limit must lie in 1..20")
         return
      end if
      if (self%batch_size < 1) then
         call fatal_error(error, "Projection batch size must be positive")
         return
      end if
      valid = ieee_is_finite(self%solve%switch_fnorm)
      if (valid) valid = self%solve%switch_fnorm >= 0.0_wp
      if (.not. valid) then
         call fatal_error(error, "Newton stage-switch residual must be nonnegative")
         return
      end if
      if (self%solve%max_cheap_steps < 0) then
         call fatal_error(error, "Newton cheap-step budget must not be negative")
         return
      end if
      valid = ieee_is_finite(self%solve%max_step)
      if (valid) valid = self%solve%max_step > 0.0_wp
      if (.not. valid) then
         call fatal_error(error, "Newton step cap must be positive")
         return
      end if
      valid = ieee_is_finite(self%solve%merit_margin)
      if (valid) valid = self%solve%merit_margin > 0.0_wp
      if (.not. valid) then
         call fatal_error(error, "Newton merit margin must be positive")
         return
      end if

      !> Derive tolerance hierarchy from master tolerance
      !> (tightest to loosest: wleb_cut < screening < solve tol < branch_sep)
      self%finalize%wleb_cut = self%tolerance*0.05_wp
      self%screening_threshold = self%tolerance*0.1_wp
      self%solve%tol = self%tolerance
      ! Copies of an ill-conditioned projected point differ by ~1e-9 Bohr, distinct ones by far more
      self%branch%sep_cut = max(self%tolerance*10.0_wp, self%branch%sep_floor)

      !> Derive weight switching bounds from level
      select case (self%wleb_prune_level)
      case (0)
         self%finalize%wleb_prune_off = 0.0_wp
         self%finalize%wleb_prune_on = 0.0_wp
      case (1)
         self%finalize%wleb_prune_off = 1.0E-12_wp
         self%finalize%wleb_prune_on = 1.0E-10_wp
      case (2)
         self%finalize%wleb_prune_off = 1.0E-10_wp
         self%finalize%wleb_prune_on = 1.0E-8_wp
      case (3)
         self%finalize%wleb_prune_off = 1.0E-8_wp
         self%finalize%wleb_prune_on = 1.0E-6_wp
      case (4)
         self%finalize%wleb_prune_off = 1.0E-6_wp
         self%finalize%wleb_prune_on = 1.0E-4_wp
      case (5)
         self%finalize%wleb_prune_off = 1.0E-4_wp
         self%finalize%wleb_prune_on = 1.0E-2_wp
      case (6)
         self%finalize%wleb_prune_off = 1.0E-2_wp
         self%finalize%wleb_prune_on = 1.0E0_wp
      case default
         call fatal_error(error, "Invalid wleb_prune_level (must be 0-6)")
         return
      end select
      self%finalize%wleb_prune = self%wleb_prune_level > 0

      ! The cut stays at the tolerance with pruning on: the prune switch is
      ! exactly zero below its onset, and a raised cut would step

      !> Branch window
      !>
      !> A negative window or a non-positive alpha makes the consumers of the
      !> admissible radius take the square root of a negative number
      valid = ieee_is_finite(self%phi_alpha)
      if (valid) valid = self%phi_alpha > 0.0_wp
      if (.not. valid) then
         call fatal_error(error, "Objective alpha must be finite and positive")
         return
      end if
      valid = ieee_is_finite(self%branch%dpdist_max)
      if (valid) valid = self%branch%dpdist_max >= 0.0_wp
      if (.not. valid) then
         call fatal_error(error, "Branch distance window must not be negative")
         return
      end if
      !> The closed edge bounds the search, so it has no off switch
      if (.not. (ieee_is_finite(self%branch%reach_off) .and. self%branch%reach_on > 0.0_wp &
                 .and. self%branch%reach_off > self%branch%reach_on)) then
         call fatal_error(error, "Reach bounds must satisfy 0 < on < off")
         return
      end if
      !> Switch bands rise from off to on with nonzero width: reversed edges
      !> invert a gate, equal ones make it a step
      if (.not. (self%w_0ls_off >= 0.0_wp .and. self%w_0ls_on > self%w_0ls_off &
                 .and. ieee_is_finite(self%w_0ls_on))) then
         call fatal_error(error, "Critical-gradient switch bounds must satisfy 0 <= off < on")
         return
      end if
      !> A negative off-edge would open the gate past the focal point, det B = 0
      if (.not. (self%w_0tra_off >= 0.0_wp .and. self%w_0tra_on > self%w_0tra_off &
                 .and. ieee_is_finite(self%w_0tra_on))) then
         call fatal_error(error, "Focal switch bounds must satisfy 0 <= off < on")
         return
      end if
      if (.not. (ieee_is_finite(self%w_0ls_p) .and. ieee_is_finite(self%w_0ls_a) &
                 .and. self%w_0ls_p > 0.0_wp .and. self%w_0ls_a > 0.0_wp)) then
         call fatal_error(error, "Switch shape parameters w_0ls_p and w_0ls_a must be positive")
         return
      end if
      !> A zero on-edge switches the horizon factor off, whatever the off-edge
      if (.not. (self%w_0hor_on == 0.0_wp .or. (self%w_0hor_off >= 0.0_wp &
                                                .and. self%w_0hor_on > self%w_0hor_off &
                                                .and. self%w_0hor_on <= 1.0_wp))) then
         call fatal_error(error, "Horizon switch bounds must satisfy 0 <= off < on <= 1")
         return
      end if
      if (self%w_0hor_on > 0.0_wp) then
         self%hor_off = self%w_0hor_off
         self%hor_on = self%w_0hor_on
      else
         ! Switched off: a band below every cosine in [-1, 1], factor one, slope zero
         self%hor_off = -2.0_wp
         self%hor_on = -1.0_wp
      end if

      !> Grid point adjacency list cutoff from density kernel length
      self%adj_list_grid_cutoff = 4.0_wp*self%rho_grid_h

      !> `do_fine` switches every optional property on; never switches one off
      if (self%do_fine) then
         self%do_grid_density = .true.
         self%do_curvature = .true.
         self%do_normal = .true.
         self%do_atomarea = .true.
      end if

      !> Select fitted Born zeta for the active Lebedev grid size
      call self%select_born_zeta(error)

   end subroutine compute_drop_derived

   !> Register all user-facing parameter entries for JSON/TOML configuration
   !>
   !> Connects parameter fields to their dotted key names for
   !> automatic parsing from configuration files via read_file/write_file
   !> Derived parameters are not registered (they are recomputed)
   !>
   !> @param[inout] self Parameter container
   subroutine register_cavity_drop_entries(self)
      class(moist_cavity_drop_parameters_type), intent(inout), target :: self

      call self%register_logical("do_fine", self%do_fine)
      ! Optional properties
      call self%register_logical("do_grid_density", self%do_grid_density)
      call self%register_logical("do_curvature", self%do_curvature)
      call self%register_logical("do_normal", self%do_normal)
      call self%register_logical("do_atomarea", self%do_atomarea)
      ! Grid
      call self%register_int_scalar("grid.num_leb", self%num_leb)
      ! Tolerance (master; wleb_cut, solve tol, branch_sep are derived;
      ! LSF screening lives on the LSF concrete and is not registered here)
      call self%register_real_scalar("tolerance", self%tolerance)
      ! Objective
      call self%register_real_scalar("objective.alpha", self%phi_alpha)
      ! Projection
      call self%register_int_scalar("projection.maxiter", self%solve%max_iter)
      call self%register_int_scalar("projection.level", self%proj_level)
      ! Octree surface coverage (projection levels 3 and 4, fallback of 1)
      call self%register_int_scalar("projection.octree.max_boxes", self%seed%max_boxes)
      call self%register_int_scalar("projection.octree.max_survivors", self%seed%max_survivors)
      call self%register_int_scalar("projection.octree.max_depth", self%seed%max_depth)
      ! Batched projection
      call self%register_int_scalar("projection.batch_size", self%batch_size)
      call self%register_real_scalar("projection.newton.switch_fnorm", self%solve%switch_fnorm)
      call self%register_int_scalar("projection.newton.max_cheap_steps", self%solve%max_cheap_steps)
      call self%register_real_scalar("projection.newton.max_step", self%solve%max_step)
      call self%register_real_scalar("projection.newton.merit_margin", self%solve%merit_margin)
      ! Screening
      call self%register_int_scalar("screening.cell_grid_full_scan_below", &
                                    self%cell_grid_full_scan_below)
      ! Switching
      call self%register_real_scalar("switching.w_0ls_off", self%w_0ls_off)
      call self%register_real_scalar("switching.w_0ls_on", self%w_0ls_on)
      call self%register_real_scalar("switching.w_0ls_p", self%w_0ls_p)
      call self%register_real_scalar("switching.w_0ls_a", self%w_0ls_a)
      call self%register_real_scalar("switching.w_0tra_off", self%w_0tra_off)
      call self%register_real_scalar("switching.w_0tra_on", self%w_0tra_on)
      call self%register_real_scalar("switching.w_0hor_off", self%w_0hor_off)
      call self%register_real_scalar("switching.w_0hor_on", self%w_0hor_on)
      ! Weight switching (off/on are derived)
      call self%register_int_scalar("switching.wleb_prune_level", self%wleb_prune_level)
      ! Density
      call self%register_real_scalar("density.rho_grid_h", self%rho_grid_h)
      ! Branching
      call self%register_real_scalar("branching.dpdist_max", self%branch%dpdist_max)
      call self%register_real_scalar("branching.reach_on", self%branch%reach_on)
      call self%register_real_scalar("branching.reach_off", self%branch%reach_off)

   end subroutine register_cavity_drop_entries

   !> Load parameters from a JSON or TOML file and recompute derived values
   !>
   !> Wraps the inherited read_file (ensure_entries -> init_defaults ->
   !> read a file), then recomputes derived parameters
   !>
   !> @param[inout] self Parameter container
   !> @param[in]    filepath Path to the JSON or TOML parameter file
   !> @param[out]   error Error if no fitted Born zeta exists for loaded num_leb
   subroutine load_drop_file(self, filepath, error)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      character(len=*), intent(in) :: filepath
      type(error_type), allocatable, intent(out) :: error

      call self%read_file(filepath, error)
   end subroutine load_drop_file

   !> Select the fitted iSwiG Born zeta parameter for the active Lebedev grid
   !>
   !> @param[inout] self Parameter container
   !> @param[out]   error Error if no fitted value exists for `num_leb`
   subroutine select_born_zeta(self, error)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      logical :: found

      call get_born_zeta(self%num_leb, self%iswig_xi_born, found)
      if (.not. found) call fatal_error(error, "No fitted DROP Born zeta parameter for Lebedev grid")
   end subroutine select_born_zeta

   !> Whether `level` is a supported projection level: 1, 3 or 4
   !>
   !> @param[in] level Projection level
   !> @return    valid Supported level
   pure function is_valid_proj_level(level) result(valid)
      integer, intent(in) :: level
      logical :: valid

      valid = level == 1 .or. level == 3 .or. level == 4
   end function is_valid_proj_level

   !> Return a human-readable label for the current projection level
   !>
   !> @param[in] self Parameter container
   pure function get_proj_level_label(self) result(label)
      !> Parameter container
      class(moist_cavity_drop_parameters_type), intent(in) :: self
      !> Description of the active projection level, blank-padded
      character(len=proj_level_label_len) :: label

      select case (self%proj_level)
      case (1)
         label = "Local Newton, octree retry"
      case (3)
         label = "Coarse octree seeds"
      case (4)
         label = "Fine octree seeds"
      case default
         label = "unknown"
      end select
   end function get_proj_level_label

   !> Print current parameter values
   !>
   !> @param[in] self Parameter container to display
   !> @param[in] unit Output unit (default `output_unit`); callers holding a run
   !>                 context pass `ctx%unit` so this honours a log file
   !> @param[in] level Verbosity level; below 3 only the normal settings print
   !>                  (default: all)
   subroutine print_parameters(self, unit, level)
      class(moist_cavity_drop_parameters_type), intent(in) :: self
      !> Output unit override
      integer, intent(in), optional :: unit
      !> Verbosity level of the printout
      integer, intent(in), optional :: level
      type(prettyprinter) :: pp
      !> Effective output unit
      integer :: iu
      !> Whether the detailed settings print
      logical :: detailed

      iu = output_unit
      if (present(unit)) iu = unit
      detailed = .true.
      if (present(level)) detailed = level >= 3

      pp = new_prettyprinter(unit=iu)

      call pp%push("DROP Parameters:")

      call pp%push("Discretization:")
      call pp%kv("Number of Leb. points", self%num_leb)
      call pp%pop()

      call pp%push("Tolerance:")
      call pp%kv("Main tolerance", self%tolerance)
      if (detailed) then
         call pp%kv("Projection tol.", self%solve%tol)
         call pp%kv("Screening threshold", self%screening_threshold)
         call pp%kv("Weight cutoff", self%finalize%wleb_cut)
         call pp%kv("Branch sep. cutoff", self%branch%sep_cut)
      end if
      call pp%pop()

      if (detailed) then
         call pp%push("Switching:")
         call pp%kv("f_crit off", self%w_0ls_off)
         call pp%kv("f_crit on", self%w_0ls_on)
         call pp%kv("f_foc off", self%w_0tra_off)
         call pp%kv("f_foc on", self%w_0tra_on)
         call pp%kv("f_hor off", self%w_0hor_off)
         call pp%kv("f_hor on", self%w_0hor_on)
         call pp%kv("Wleb switch level", self%wleb_prune_level)
         if (self%finalize%wleb_prune) then
            call pp%kv("Wleb switch off", self%finalize%wleb_prune_off)
            call pp%kv("Wleb switch on", self%finalize%wleb_prune_on)
         end if
         call pp%pop()

         call pp%push("Gaussians:")
         call pp%kv("xi_born", self%iswig_xi_born)
         call pp%pop()

         call pp%push("Objective function:")
         call pp%kv("Alpha", self%phi_alpha)
         call pp%pop()
      end if

      call pp%push("Projection settings:")
      call pp%kv("Projection level", self%proj_level, self%proj_level_label())
      if (detailed) then
         call pp%kv("Maximum iterations", self%solve%max_iter)
         call pp%kv("Octree seed size", self%seed%seed_size, "Bohr")
         if (self%proj_level >= 3) then
            call pp%kv("Octree box budget", self%seed%max_boxes)
         end if
         call pp%kv("Batch size", self%batch_size)
      end if
      call pp%pop()

      if (detailed) then
         call pp%push("Branching:")
         call pp%kv("Max. distance excess", self%branch%dpdist_max, "Bohr")
         call pp%kv("Reach on", self%branch%reach_on, "Bohr")
         call pp%kv("Reach off", self%branch%reach_off, "Bohr")
         call pp%pop()

         call pp%push("Screening:")
         call pp%kv("Cell grid full-scan below", self%cell_grid_full_scan_below, "atoms")
         call pp%pop()

         call pp%push("Optional properties:")
         call pp%kv("Grid point density", self%do_grid_density)
         call pp%kv("Curvature", self%do_curvature)
         call pp%kv("Surface normals", self%do_normal)
         call pp%kv("Atom areas", self%do_atomarea)
         call pp%pop()
      end if

   end subroutine print_parameters

end module moist_cavity_drop_parameters
