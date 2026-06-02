module moist_cavity_drop_parameters
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use moist_model_parameters, only: moist_model_parameters_type
   use, intrinsic :: iso_fortran_env, only: output_unit
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter

   implicit none

   !> Maximum supported atomic number
   integer, parameter :: maxAtomicNumbers = 118

   public :: moist_cavity_drop_parameters_type
   private

   !> Lebedev grid sizes with available Born zeta values
   integer, parameter :: iswig_xi_born_nleb(29) = [ &
                         6, 14, 26, 38, 50, 86, 110, 146, &
                         170, 194, 302, 350, 434, 590, 770, 974, &
                         1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, &
                         3890, 4334, 4802, 5294, 5810]

   !> Fitted Born zeta values matching `iswig_xi_born_nleb`.
   real(wp), parameter :: iswig_xi_born_zeta(29) = [ &
                          4.845184_wp, 4.864049_wp, 4.854249_wp, 4.900523_wp, &
                          4.891966_wp, 4.896867_wp, 4.900490_wp, 4.897689_wp, &
                          4.906299_wp, 4.902816_wp, 4.904427_wp, 4.868267_wp, &
                          4.905109_wp, 4.905704_wp, 4.906018_wp, 4.906299_wp, &
                          4.906493_wp, 4.906667_wp, 4.906807_wp, 4.906894_wp, &
                          4.906981_wp, 4.907088_wp, 4.907121_wp, 4.907175_wp, &
                          4.907208_wp, 4.907262_wp, 4.906948_wp, 4.907088_wp, &
                          4.907402_wp]

   !> Parameter container for DROP
   !> Extends the base parameter type with DROP-specific settings for
   !> grid generation, blending, barrier potentials, and optimization.
   type, extends(moist_model_parameters_type) :: moist_cavity_drop_parameters_type

      !> ========== Grid Discretization ==========

      !> Number of Lebedev quadrature points per atomic sphere
      integer :: num_leb = 194

      !> ========== Tolerance ==========

      !> Main tolerance controlling all numerical thresholds.
      !>
      !> Derived tolerances (tightest to loosest):
      !>  - wleb_cut            = tolerance * 0.05 (quadrature weight cutoff)
      !>  - screening_threshold = tolerance * 0.1  (LSF screening; passed into lsf_model by DROP constructor)
      !>  - proj_tol            = tolerance        (projection convergence)
      !>  - branch_sep_cut      = tolerance * 10   (branch degeneracy)
      real(wp) :: tolerance = 1.0E-10_wp

      !> Minimum weight cutoff (derived from tolerance)
      real(wp) :: wleb_cut = 5.0E-12_wp

      !> LSF screening threshold (derived from `tolerance` in
      !> `compute_drop_derived`). The DROP constructor pushes this value
      !> into the polymorphic `lsf_model` so the LSF's own internal
      !> screening caches use a value consistent with the cavity tolerance.
      real(wp) :: screening_threshold = 1.0E-11_wp

      !> ========== Objective Function Weights ==========

      !> Weight w_a for anchor term (keeps points near initial position)
      real(wp) :: phi_alpha = 0.5_wp

      !> ========== Projection Optimization ==========

      !> Convergence tolerance for projection optimizer (derived from tolerance)
      real(wp) :: proj_tol = 1.0E-10_wp
      !> Maximum iterations for projection optimizer
      integer :: proj_maxiter = 150
      !> Projection level:
      !> 1 = SLSQP only (no Newton refinement)
      !> 2 = SLSQP + Newton refinement
      !> 3 = Conditional multi-tangent for degenerate points
      !> 4 = Conditional SLSQP-deflation
      !> 5 = SLSQP-deflation (unconditional)
      !> 6 = Newton-deflation on 4D KKT system
      !> 7 = Regular SLSQP multistart
      !> 8 = Fine SLSQP multistart reference (very expensive)
      integer :: proj_level = 3

      !> ========== Screening ==========

      !> Distance cutoff for grid point adj. list
      real(wp) :: adj_list_grid_cutoff = 1.0_wp
      !> Below this atom count, the cell grid collapses to a single
      !> full-scan cell. For small systems the per-cell fan-out reduces
      !> to "every atom in every cell" anyway, so we skip the build work
      !> and let every query return the full atom list directly
      integer :: cell_grid_full_scan_below = 200
      !> Cell fraction for molecular cell grid
      !> (1.0 = no subdivision, 0.5 = halved cell,..)
      real(wp) :: cell_grid_fraction = 0.25_wp

      !> ========== Hard-sphere reference cavity ========

      !> Active iSwiG Gaussian zeta fitted to the analytical Born energy
      !> for the selected Lebedev grid size
      real(wp) :: iswig_xi_born

      !> ========== Switching Functions ==========

      !> Start of critical level-set weight switching transition
      real(wp) :: w_0ls_from = 0.25_wp
      !> End of critical level-set weight switching transition
      real(wp) :: w_0ls_to = 0.6_wp

      !> Start of focal/branching weight switching transition
      real(wp) :: w_0tra_from = 0.1_wp
      !> End of focal/branching weight switching transition
      real(wp) :: w_0tra_to = 0.3_wp

      !> Function parameters
      real(wp) :: w_0ls_p = 0.8_wp
      real(wp) :: w_0ls_a = 1.6_wp

      !> ========== Grid point density ==========

      !> Grid density kernel length
      real(wp) :: rho_grid_h = 1.0_wp

      !> ========== Branching ==========

      !> Softmax scale parameter for branch weight model (smoothness)
      real(wp) :: branch_weight_s = 0.0025_wp
      !> Rho cutoff for branch separation (derived from `branch_weight_s` and `wleb_cut`)
      real(wp) :: branch_rho_cut
      !> Branch separation cutoff (derived from tolerance)
      real(wp) :: branch_sep_cut = 1.0E-8_wp

      !> ========== Weight switching ==========

      !> Smooth switching on final Lebedev weights to suppress near-zero contributions before the branch filter.
      !>  - Level 0: disabled (default)
      !>  - Level 1: from 1E-12 to 1E-10
      !>  - Level 2: from 1E-10 to 1E-8
      !>  - Level 3: from 1E-8  to 1E-6
      !>  - Level 4: from 1E-6  to 1E-4
      !>  - Level 5: from 1E-4  to 1E-2
      !>  - Level 6: from 1E-2  to 1E0
      integer :: wleb_prune_level = 0
      !> Lower bound of the weight switching region (below: fully off, derived)
      real(wp) :: wleb_prune_from = 0.0_wp
      !> Upper bound of the weight switching region (above: fully on, derived)
      real(wp) :: wleb_prune_to = 0.0_wp

      !> ========== Disconnected points ==========

      !> Point disconnection distance threshold (times the average grid point spacing)
      real(wp) :: disconnection_thrs = 4.0_wp

   contains
      !> Initialize parameters to compiled defaults
      procedure :: init_defaults => init_cavity_drop_defaults
      !> Initialize parameters from constructor inputs.
      procedure :: new => new_moist_cavity_drop_parameters_type
      !> Register parameters for JSON configuration parsing
      procedure :: register_entries => register_cavity_drop_entries
      !> Recompute all derived parameters from user-facing fields
      procedure, private :: compute_derived => compute_drop_derived
      !> Load parameters from JSON file and recompute derived values.
      !> Use this instead of read_file to ensure derived parameters
      !> (Born zeta, Jacobian regularization, etc.) stay consistent.
      procedure :: load_file => load_drop_file
      !> Select the fitted Born zeta value for the active Lebedev grid.
      procedure :: select_born_zeta => select_born_zeta
      !> Print current parameter values to output
      procedure :: print => print_parameters
      !> Return a human-readable description for the current projection level
      procedure :: proj_level_label => get_proj_level_label

   end type moist_cavity_drop_parameters_type

contains

   !> Construct a new DROP parameters instance
   !> @param[inout] self Initialized parameter container
   !> @param[in]    nleb Number of Lebedev quadrature points
   !> @param[in]    tolerance main tolerance (derives proj_tol, wleb_cut, branch_sep)
   !> @param[in]    proj_maxiter Maximum projection iterations
   !> @param[in]    proj_level Projection refinement level
   !> @param[in]    branch_weight_s Softmax scale for branch weights
   !> @param[in]    rho_grid_h Grid density kernel length
   !> @param[in]    wleb_prune_level Weight switching level (0=off, 1-6=increasing aggressiveness)
   !> @param[out]   error Error if no fitted Born zeta exists for `nleb`
   subroutine new_moist_cavity_drop_parameters_type(self, &
                                                    nleb, tolerance, proj_maxiter, proj_level, &
                                                    branch_weight_s, rho_grid_h, wleb_prune_level, error)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      integer, intent(in), optional :: nleb
      real(wp), intent(in), optional :: tolerance
      integer, intent(in), optional :: proj_maxiter
      integer, intent(in), optional :: proj_level
      real(wp), intent(in), optional :: branch_weight_s
      real(wp), intent(in), optional :: rho_grid_h
      integer, intent(in), optional :: wleb_prune_level
      type(error_type), allocatable, intent(out) :: error

      ! Phase 1: Reset all fields to compiled defaults
      call self%init_defaults()

      ! Phase 2: Apply caller overrides
      if (present(nleb)) self%num_leb = nleb

      if (present(tolerance)) self%tolerance = tolerance
      if (present(proj_maxiter)) self%proj_maxiter = proj_maxiter
      if (present(proj_level)) self%proj_level = proj_level

      if (present(branch_weight_s)) self%branch_weight_s = branch_weight_s
      if (present(rho_grid_h)) self%rho_grid_h = rho_grid_h
      if (present(wleb_prune_level)) self%wleb_prune_level = wleb_prune_level

      ! Phase 3: Compute derived parameters from final field values
      call self%compute_derived(error)
      if (allocated(error)) return

   end subroutine new_moist_cavity_drop_parameters_type

   !> Reset all user-facing parameters to compiled default values.
   !> Defaults are defined in the type declaration; this is called by
   !> read_file before loading JSON to ensure a clean slate on re-reads.
   !> @param[inout] self Parameter container to reset
   subroutine init_cavity_drop_defaults(self)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      !> Fresh instance carrying the compiled defaults from the type declaration
      type(moist_cavity_drop_parameters_type) :: fresh

      ! Grid
      self%num_leb = fresh%num_leb
      ! Tolerance (master; proj_tol, wleb_cut, branch_sep_cut are derived
      ! from this in compute_derived; LSF screening lives on the LSF concrete)
      self%tolerance = fresh%tolerance
      ! Objective
      self%phi_alpha = fresh%phi_alpha
      ! Projection
      self%proj_maxiter = fresh%proj_maxiter
      self%proj_level = fresh%proj_level
      ! Screening
      self%cell_grid_full_scan_below = fresh%cell_grid_full_scan_below
      self%cell_grid_fraction = fresh%cell_grid_fraction
      ! Switching
      self%w_0ls_from = fresh%w_0ls_from
      self%w_0ls_to = fresh%w_0ls_to
      self%w_0ls_p = fresh%w_0ls_p
      self%w_0ls_a = fresh%w_0ls_a
      self%w_0tra_from = fresh%w_0tra_from
      self%w_0tra_to = fresh%w_0tra_to
      ! Weight switching (from/to are derived in compute_derived)
      self%wleb_prune_level = fresh%wleb_prune_level
      ! Density
      self%rho_grid_h = fresh%rho_grid_h
      ! Branching
      self%branch_weight_s = fresh%branch_weight_s
      ! Disconnected points
      self%disconnection_thrs = fresh%disconnection_thrs
   end subroutine init_cavity_drop_defaults

   !> Compute all derived parameters from the current user-facing fields
   !>
   !> Must be called after any modification of user-facing parameters
   !> (constructor or file load) to keep derived values consistent
   !> @param[inout] self Parameter container
   !> @param[out]   error Error if no fitted Born zeta exists for `num_leb`
   subroutine compute_drop_derived(self, error)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      !> Derive tolerance hierarchy from master tolerance
      !> (tightest to loosest: wleb_cut < screening < proj_tol < branch_sep)
      self%wleb_cut = self%tolerance*0.05_wp
      self%screening_threshold = self%tolerance*0.1_wp
      self%proj_tol = self%tolerance
      self%branch_sep_cut = self%tolerance*10.0_wp

      !> Branch rho cutoff from weight cutoff and softmax scale
      self%branch_rho_cut = log(1.0_wp/self%wleb_cut)*sqrt(self%branch_weight_s)

      !> Grid point adjacency list cutoff from density kernel length
      self%adj_list_grid_cutoff = 4.0_wp*self%rho_grid_h

      !> Derive weight switching bounds from level
      select case (self%wleb_prune_level)
      case (0)
         self%wleb_prune_from = 0.0_wp
         self%wleb_prune_to = 0.0_wp
      case (1)
         self%wleb_prune_from = 1.0E-12_wp
         self%wleb_prune_to = 1.0E-10_wp
      case (2)
         self%wleb_prune_from = 1.0E-10_wp
         self%wleb_prune_to = 1.0E-8_wp
      case (3)
         self%wleb_prune_from = 1.0E-8_wp
         self%wleb_prune_to = 1.0E-6_wp
      case (4)
         self%wleb_prune_from = 1.0E-6_wp
         self%wleb_prune_to = 1.0E-4_wp
      case (5)
         self%wleb_prune_from = 1.0E-4_wp
         self%wleb_prune_to = 1.0E-2_wp
      case (6)
         self%wleb_prune_from = 1.0E-2_wp
         self%wleb_prune_to = 1.0E0_wp
      case default
         call fatal_error(error, "Invalid wleb_prune_level (must be 0-6)")
         return
      end select

      ! If prune level is on, adjust wleb_cut
      if (self%wleb_prune_level > 0) then
         self%wleb_cut = self%wleb_prune_from
      end if

      !> Select fitted Born zeta for the active Lebedev grid size
      call self%select_born_zeta(error)

   end subroutine compute_drop_derived

   !> Register all user-facing parameter entries for JSON configuration.
   !>
   !> Connects parameter fields to their dotted key names for
   !> automatic parsing from configuration files via read_file/write_file.
   !> Derived parameters are not registered (they are recomputed).
   !> @param[inout] self Parameter container
   subroutine register_cavity_drop_entries(self)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self

      ! Grid
      call self%register_int_scalar('grid.num_leb', self%num_leb)
      ! Tolerance (master; wleb_cut, proj_tol, branch_sep are derived;
      ! LSF screening lives on the LSF concrete and is not registered here)
      call self%register_real_scalar('tolerance', self%tolerance)
      ! Objective
      call self%register_real_scalar('objective.alpha', self%phi_alpha)
      ! Projection
      call self%register_int_scalar('projection.maxiter', self%proj_maxiter)
      call self%register_int_scalar('projection.level', self%proj_level)
      ! Screening
      call self%register_int_scalar('screening.cell_grid_full_scan_below', &
                                    self%cell_grid_full_scan_below)
      call self%register_real_scalar('screening.cell_grid_fraction', &
                                     self%cell_grid_fraction)
      ! Switching
      call self%register_real_scalar('switching.w_0ls_from', self%w_0ls_from)
      call self%register_real_scalar('switching.w_0ls_to', self%w_0ls_to)
      call self%register_real_scalar('switching.w_0ls_p', self%w_0ls_p)
      call self%register_real_scalar('switching.w_0ls_a', self%w_0ls_a)
      call self%register_real_scalar('switching.w_0tra_from', self%w_0tra_from)
      call self%register_real_scalar('switching.w_0tra_to', self%w_0tra_to)
      ! Weight switching (from/to are derived)
      call self%register_int_scalar('switching.wleb_prune_level', self%wleb_prune_level)
      ! Density
      call self%register_real_scalar('density.rho_grid_h', self%rho_grid_h)
      ! Branching
      call self%register_real_scalar('branching.softmax_scale', self%branch_weight_s)
      ! Disconnected points
      call self%register_real_scalar('disconnection.threshold', self%disconnection_thrs)

   end subroutine register_cavity_drop_entries

   !> Load parameters from a JSON file and recompute derived values
   !>
   !> Wraps the inherited read_file (ensure_entries -> init_defaults ->
   !> read JSON), then recomputes derived parameters
   !> @param[inout] self Parameter container
   !> @param[in]    filepath Path to the JSON parameter file
   !> @param[out]   error Error if no fitted Born zeta exists for loaded num_leb
   subroutine load_drop_file(self, filepath, error)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      character(len=*), intent(in) :: filepath
      type(error_type), allocatable, intent(out) :: error

      ! Delegate to inherited read_file: ensure_entries -> init_defaults -> read JSON
      call self%read_file(filepath)

      ! Recompute derived parameters from loaded values
      call self%compute_derived(error)
   end subroutine load_drop_file

   !> Select the fitted iSwiG Born zeta parameter for the active Lebedev grid
   !> @param[inout] self Parameter container
   !> @param[out]   error Error if no fitted value exists for `num_leb`
   subroutine select_born_zeta(self, error)
      class(moist_cavity_drop_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      integer :: igrid

      do igrid = 1, size(iswig_xi_born_nleb)
         if (self%num_leb == iswig_xi_born_nleb(igrid)) then
            self%iswig_xi_born = iswig_xi_born_zeta(igrid)
            return
         end if
      end do

      call fatal_error(error, "No fitted DROP Born zeta parameter for Lebedev grid")
   end subroutine select_born_zeta

   !> Return a human-readable label for the current projection level
   !> @param[in] self Parameter container
   !> @return    label Description of the active projection level
   pure function get_proj_level_label(self) result(label)
      class(moist_cavity_drop_parameters_type), intent(in) :: self
      character(len=:), allocatable :: label

      select case (self%proj_level)
      case (1)
         label = 'SLSQP'
      case (2)
         label = 'SLSQP + Newton'
      case (3)
         label = 'Conditional multi-tangent'
      case (4)
         label = 'Conditional SLSQP-deflation'
      case (5)
         label = 'SLSQP-deflation'
      case (6)
         label = 'Newton-deflation (4D KKT)'
      case (7)
         label = 'Regular SLSQP multistart'
      case (8)
         label = 'Fine SLSQP multistart reference'
      case default
         label = 'unknown'
      end select
   end function get_proj_level_label

   !> Print current parameter values to standard output
   !> @param[in] self Parameter container to display
   subroutine print_parameters(self)
      class(moist_cavity_drop_parameters_type), intent(in) :: self
      type(prettyprinter) :: pp

      pp = new_prettyprinter(unit=output_unit)

      call pp%blank()
      call pp%push('Cavity Parameters:')

      call pp%push('Discretization:')
      call pp%kv('Number of Leb. points', self%num_leb)
      call pp%pop()

      call pp%push('Tolerance:')
      call pp%kv('Main tolerance', self%tolerance)
      call pp%kv('Projection tol.', self%proj_tol)
      call pp%kv('Screening threshold', self%screening_threshold)
      call pp%kv('Weight cutoff', self%wleb_cut)
      call pp%kv('Branch sep. cutoff', self%branch_sep_cut)
      call pp%pop()

      call pp%push('Switching:')
      call pp%kv('f_crit start', self%w_0ls_from)
      call pp%kv('f_crit end', self%w_0ls_to)
      call pp%kv('f_foc start', self%w_0tra_from)
      call pp%kv('f_foc end', self%w_0tra_to)
      call pp%kv('Wleb switch level', self%wleb_prune_level)
      if (self%wleb_prune_level > 0) then
         call pp%kv('Wleb switch from', self%wleb_prune_from)
         call pp%kv('Wleb switch to', self%wleb_prune_to)
      end if
      call pp%pop()

      call pp%push('Gaussians:')
      call pp%kv('xi_born', self%iswig_xi_born)
      call pp%pop()

      call pp%push('Objective function:')
      call pp%kv('Alpha', self%phi_alpha)
      call pp%pop()

      call pp%push('Projection settings:')
      call pp%kv('Projection level', self%proj_level, self%proj_level_label())
      call pp%kv('Maximum iterations', self%proj_maxiter)
      call pp%pop()

      call pp%push('Branching:')
      call pp%kv('Softmax scale', self%branch_weight_s)
      call pp%kv('Rho cutoff', self%branch_rho_cut, 'Bohr')
      call pp%pop()

      call pp%push('Screening:')
      call pp%kv('Cell grid full-scan below', self%cell_grid_full_scan_below, 'atoms')
      call pp%kv('Cell grid fraction', self%cell_grid_fraction)
      call pp%pop()

      call pp%push('Disconnected points:')
      call pp%kv('Distance threshold', self%disconnection_thrs)
      call pp%pop()

   end subroutine print_parameters

end module moist_cavity_drop_parameters
