!> Main DROP (Discretization via Reference-Onto surface Projection) implementation
module moist_cavity_drop
   use mctc_env, only: wp
   use mctc_io_structure, only: structure_type
   use mctc_io, only: new
   use mctc_env, only: error_type, fatal_error, wp
   use moist_math_lapack, only: getrf, getrs
   use moist_math_lapack_gesv, only: dgesv
   use moist_math_linalg, only: mat3x3_inv, setup_tangent_frame
   use moist_math_boys, only: dboysfun1
   use moist_math_grid_lebedev, only: get_angular_grid, grid_size, lebedev_order_from_num
   use moist_type, only: cavity_type, potential_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_context, only: moist_context_type
   use moist_radius_type, only: radius_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_switching, only: moist_cavity_drop_smooth_step_swif, new_smooth_step_swif
   use moist_cavity_drop_switching, only: moist_cavity_drop_sigmoid_bump_swif, new_sigmoid_bump_swif
   use moist_cavity_drop_gaussian, only: moist_cavity_drop_iswig, new_iswig
   use moist_cavity_drop_projector, only: drop_projector_type
   use moist_cavity_drop_types, only: projection_buffer_type, projection_workspace_type
   use moist_math_adjacency_list, only: adjacency_list_type
   use moist_math_cell_grid, only: moist_cell_grid_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_utils_mem, only: grow_array, filter_array
   use moist_cavity_drop_objective_phi, only: moist_cavity_drop_objective_phi_type
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, drop_surface_weights_type
   use moist_math_smoothing_kernels, only: wendland_kernel_type

   use moist_utils_timer, only: timer_type, cat_setup, cat_solve, cat_properties, cat_gradient
   use moist_cavity_drop_request, only: drop_property_request, &
                                        drop_request_default, drop_request_diagnostics, drop_request_fine

   implicit none
   private

   public :: cavity_type_drop
   public :: new_cavity_drop
   public :: drop_property_request
   public :: drop_request_default, drop_request_diagnostics, drop_request_fine

   !> DROP cavity type
   type, extends(cavity_type) :: cavity_type_drop

      !> DROP parameters (blend_k, blend_1b, blend_2b, blend_3b, etc.)
      type(moist_cavity_drop_parameters_type) :: param

      !> Level set function model
      class(moist_cavity_drop_lsf_type), allocatable :: lsf_model

      !> Property request flags controlling which quantities are computed
      type(drop_property_request) :: request

      !* ----------------------------- Atomic sphere data ----------------------------- *!

      !> Molecular structure
      type(structure_type) :: mol

      !> Per-cell candidate atom lists for point-to-atom screening
      type(moist_cell_grid_type) :: mol_cell_grid

      !> Grid-point neighbour list for density computation (CSR format)
      type(adjacency_list_type) :: grid_adj_list

      !> Unique numbering
      integer, allocatable :: numbering(:)

      !* ------------------------------ Integration grid ------------------------------ *!

      !> Anchor Lebedev quadrature weights (ngrid)
      real(wp), allocatable :: anchor_wleb0(:)
      !> Final Lebedev quadrature weights (ngrid)
      real(wp), allocatable :: wleb(:)
      !> Raw Lebedev quadrature weights (3, nsph, ngrid)
      real(wp), allocatable :: wleb1_rA(:, :, :)

      !> Grid-point density scaling
      real(wp), allocatable :: rho_scal0(:)

      !> Closest point Jacobian
      real(wp), allocatable :: cpjac_scal0(:)
      !> Per-gridpoint derivatives of cpjac_scal w.r.t. nuclear coordinates (3, nsph, ngrid)
      real(wp), allocatable :: cpjac_scal1_rA(:, :, :)

      !> Local grid point density
      real(wp), allocatable :: rho_grid(:)
      !> Local grid point density (hard sphere limit)
      real(wp), allocatable :: rho_grid_anchor(:)

      !* ---------------------------- Geometric quantities ---------------------------- *!

      !> Initial (unprojected) grid positions
      real(wp), allocatable :: anchorxyz(:, :)

      !> Per-gridpoint second derivatives of gridpoint w.r.t. nuclear coordinates
      !> Shape: (3, 3, nsph, 3, nsph, ngrid) = (j, alpha, A, beta, B, igrid)
      real(wp), allocatable :: xyz2_rArB(:, :, :, :, :, :)

      !> Distance from anchor to projected point
      real(wp), allocatable :: rho(:)
      !> Per-gridpoint derivatives of rho w.r.t. nuclear coordinates (3, nsph, ngrid)
      real(wp), allocatable :: rho1_rA(:, :, :)

      !> Distance from sphere center to grid point
      real(wp), allocatable :: r_iI0(:)
      !> Gridpoint-owner distance derivatives r_iI1_rA (3, nsph, ngrid)
      real(wp), allocatable :: r_iI1_rA(:, :, :)

      !> Surface normal gradient at grid points (3, nsph, 3, ngrid)
      real(wp), allocatable :: normal1_rA(:, :, :, :)

      !> First principal curvature
      real(wp), allocatable :: k1(:)
      !> Second principal curvature
      real(wp), allocatable :: k2(:)
      !> Mean curvature
      real(wp), allocatable :: KM(:)
      !> Gaussian curvature
      real(wp), allocatable :: KG(:)

      !> Nuclear gradient of the first principal curvature (3, nsph, ngrid).
      !> Mean/Gaussian curvature gradients are derived from k1_rA/k2_rA
      !> downstream: dKM = (k1_rA + k2_rA)/2, dKG = k2*k1_rA + k1*k2_rA.
      real(wp), allocatable :: k1_rA(:, :, :)
      !> Nuclear gradient of the second principal curvature (3, nsph, ngrid)
      real(wp), allocatable :: k2_rA(:, :, :)

      !* ------------------------------------ Areas ----------------------------------- *!

      !> Per-gridpoint derivatives of the area (3, nsph, ngrid)
      real(wp), allocatable :: a_i1_rA(:, :, :)
      !> Accumulated gradient of per-atom areas w.r.t. nuclear coordinates (3, nsph, nsph_owner)
      real(wp), allocatable :: asph1_rA(:, :, :)
      !> Accumulated gradient of total area w.r.t. nuclear coordinates (3, nsph)
      real(wp), allocatable :: A_tot1_rA(:, :)

      !* ----------------------------------- Volumes ---------------------------------- *!

      !> Atomic volumes
      real(wp), allocatable :: vsph(:)
      !> Accumulated gradient of per-atom volumes w.r.t. nuclear coordinates (3, nsph, nsph)
      real(wp), allocatable :: vsph1_rA(:, :, :)
      !> Accumulated gradient of total volume w.r.t. nuclear coordinates (3, nsph)
      real(wp), allocatable :: V_tot1_rA(:, :)

      !* -------------------------- Grid settings and caching ------------------------- *!

      !> Total number of initial grid points (before filtering)
      integer, allocatable :: nmax
      !> Lebedev grid order index
      integer, allocatable :: oleb
      !> Cached Lebedev angular grid (3, num_leb)
      real(wp), allocatable :: ang_grid(:, :)
      !> Cached Lebedev weights (num_leb)
      real(wp), allocatable :: ang_weight(:)

      !* ----------------------------- Switching functions ---------------------------- *!

      !> Critical level set switching function
      type(moist_cavity_drop_sigmoid_bump_swif) :: f_crit
      !> Focal/branching point switching function
      type(moist_cavity_drop_sigmoid_bump_swif) :: f_foc
      !> Lebedev weight pruning function (suppresses near-zero weights)
      type(moist_cavity_drop_sigmoid_bump_swif) :: f_wleb

      !> Weight switching values per point (ngrid)
      real(wp), allocatable :: w_f0(:)
      !> Individual switching function component derivatives (3, nsph, ngrid)
      real(wp), allocatable :: w_f1_rA(:, :, :)

      !> iSwiG
      type(moist_cavity_drop_iswig) :: iswig

      !> iSwiG switching values per point
      real(wp), allocatable :: iswig_f0(:)

      !* --------------------------------- Projection --------------------------------- *!

      !> Branch weight model for handling degenerate branches
      type(branch_weight_type) :: branch_weight

      !> Branch index per grid point (1 = default/unbranched)
      integer, allocatable :: branch(:)
      !> Anchor group id for each grid point
      integer, allocatable :: anchor_id(:)
      !> Number of branches in the anchor group of each grid point
      integer, allocatable :: branch_count(:)
      !> Per-grid-point projection convergence flag (ngrid)
      logical, allocatable :: converged(:)
      !> Branch weight per grid point (ngrid) (1 = default/unbranched)
      real(wp), allocatable :: wbranch(:)
      !> Projection-objective value phi at each surviving grid point (ngrid).
      !> Used by the branch-weight softmax gradient to recover per-branch
      !> phi values once the main-loop KKT solve has produced dr/dR_A.
      real(wp), allocatable :: phi0(:)

      real(wp), allocatable :: lambda0(:)
      !> Per-gridpoint derivatives of lambda w.r.t. nuclear coordinates (3, nsph, ngrid)
      real(wp), allocatable :: lambda1_rA(:, :, :)
      !> Per-gridpoint second derivatives of lambda w.r.t. nuclear coordinates
      !> Shape: (3, nsph, 3, nsph, ngrid) = (alpha, A, beta, B, igrid)
      real(wp), allocatable :: lambda2_rArB(:, :, :, :, :)

      !* ------------------------------------ CPCM ------------------------------------ *!

      !> Gaussian width at each anchor grid point (ngrid)
      real(wp), allocatable :: anchor_xi0(:)

   contains
      !> Configure which optional properties to compute
      procedure :: properties => set_properties_drop
      !> Update cavity for new geometry
      procedure :: update => update_cavity_drop
      !> Compute area gradient w.r.t. nuclear coordinates
      procedure :: get_gradient => get_gradient_drop

      !> Compute CPCM properties
      procedure :: compute_gaussians

      !> Contract surface-coordinate weights to LSF adjoint weights
      procedure :: contract_surface_lsf_weights
      !> Map surface-coordinate weights into generic potential channels
      procedure :: get_surface_potential => get_surface_potential_drop
      !> Contract surface-coordinate weights into the nuclear gradient
      procedure :: get_surface_gradient => get_surface_gradient_drop

      !> Compute all needed cavity gradients
      procedure :: compute_gradient_drop

      !> Compute anchor-only nuclear derivatives (callback/isodensity LSF)
      procedure :: compute_anchor_gradient

      !> Build per-cell atom lists for point-to-atom screening
      procedure :: setup_mol_cell_grid
      !> Build grid-point neighbour list for density computation
      procedure :: setup_grid_adj_list

      !> Initialize grid arrays with Lebedev points on spheres
      procedure :: fill_arrays

      !> Switching functions
      procedure :: compute_switching => compute_switching_function

      !> Remove points below switching cutoff (after projection)
      procedure :: filter_arrays

      !> Compute branch weights on the final surviving set and fold
      !> them into wleb. Runs once, after all filter passes.
      procedure :: compute_branch_weights

      !> Compute surface area from projected grid
      procedure :: compute_area_volume

      !> Cache Lebedev grid for current num_leb
      procedure :: ensure_lebedev_cache
      !> Project all grid points onto LSF surface
      procedure :: project_all_points
      !> Compute mean and Gaussian curvatures
      procedure :: compute_curvature
      !> Compute closest-point Jacobian scaling
      procedure :: compute_cp_jacobian_scaling
      !> Compute local grid point densities
      procedure :: compute_grid_point_density
      !> Run diagnostic checks on grid point health
      procedure :: analyze_cavity
      !> Write grid to CSV file for analysis
      procedure :: write_csv_debug => write_cavity_csv_debug

      !> Finalizer
      final :: finalize_cavity_drop

   end type cavity_type_drop

   interface

      !* ============================================================================== *!
      !*               Internal DROP routines (should not be used outside)               *!
      !* ============================================================================== *!

      !> [projection.f90] Compute the next capacity for projection work arrays
      module pure integer function projection_grow_capacity(current_capacity, required_capacity) result(new_capacity)
         implicit none (type, external)
         integer, intent(in) :: current_capacity
         integer, intent(in) :: required_capacity
      end function projection_grow_capacity

      !> [projection.f90] Ensure projection arrays are allocated for `new_capacity`
      module subroutine ensure_projection_capacity(self, new_capacity, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         integer, intent(in) :: new_capacity
         type(error_type), allocatable, intent(out) :: error
      end subroutine ensure_projection_capacity

      !> [projection.f90] Project all current grid points onto the LSF surface
      module subroutine project_all_points(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine project_all_points

      !> [setup.f90] Evaluate switching weights for all grid points
      module subroutine compute_switching_function(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_switching_function

      !> [properties.f90] Compute local grid-point density values
      module subroutine compute_grid_point_density(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_grid_point_density

      !> [properties.f90] Check cavity diagnostics
      module subroutine analyze_cavity(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine analyze_cavity

      !> [setup.f90] Fill initial DROP arrays from per-atom Lebedev grids
      module subroutine fill_arrays(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine fill_arrays

      !> [setup.f90] Build per-cell atom screening grid
      module subroutine setup_mol_cell_grid(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine setup_mol_cell_grid

      !> [setup.f90] Build grid-point neighbour list for density computation
      module subroutine setup_grid_adj_list(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine setup_grid_adj_list

      !> [filter.f90] Filter and compact points after projection
      module subroutine filter_arrays(self, name, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         character(len=*), intent(in) :: name
         type(error_type), allocatable, intent(out) :: error
      end subroutine filter_arrays

      !> [filter.f90] Compute branch weights on the final surviving set
      module subroutine compute_branch_weights(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_branch_weights

      !> [properties.f90] Compute total/atomic area and volume contributions
      module subroutine compute_area_volume(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_area_volume

      !> [properties.f90] Compute mean and Gaussian curvature on the surface grid
      module subroutine compute_curvature(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_curvature

      !> [projection.f90] Compute closest-point Jacobian scaling factors.
      module subroutine compute_cp_jacobian_scaling(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_cp_jacobian_scaling

      !> [projection.f90] Compute CPCM Gaussian width parameters for grid points
      module subroutine compute_gaussians(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_gaussians

      !> [deriv/forward.f90] Compute first nuclear derivatives for DROP quantities
      module subroutine compute_gradient_drop(self, error, anchor_only)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
         !> Restrict each grid point's active atom to its owner (anchor motion
         !> only); used for callback/isodensity LSFs whose field nuclear
         !> derivatives are identically zero.
         logical, intent(in), optional :: anchor_only
      end subroutine compute_gradient_drop

      !> [deriv/forward.f90] Compute anchor-only nuclear derivatives (owner motion,
      !> frozen field).  Thin wrapper over compute_gradient_drop(anchor_only=.true.).
      module subroutine compute_anchor_gradient(self, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_anchor_gradient

      !> [deriv/weights.f90] Reject a surface adjoint the cavity cannot contract
      !>
      !> @param[in]  self    DROP cavity instance
      !> @param[in]  acc     Accumulated surface-observable adjoints
      !> @param[in]  context Calling routine, used to prefix the diagnostics
      !> @param[out] error   Error object
      module subroutine check_surface_adjoint(self, acc, context, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(in) :: self
         type(cavity_surface_adjoint_type), intent(in) :: acc
         character(len=*), intent(in) :: context
         type(error_type), allocatable, intent(out) :: error
      end subroutine check_surface_adjoint

      !> [deriv/weights.f90] Fold the derived weight channels and run the branch
      !> reverse pass
      !>
      !> @param[in]  self           DROP cavity instance
      !> @param[in]  acc            Accumulated surface-observable adjoints
      !> @param[in]  fold_switching Whether to fold the area channel into `w_f`
      !> @param[out] eff            Folded weights and the branch adjoint
      module subroutine prepare_surface_weights(self, acc, fold_switching, eff)
         implicit none (type, external)
         class(cavity_type_drop), intent(in) :: self
         type(cavity_surface_adjoint_type), intent(in) :: acc
         logical, intent(in) :: fold_switching
         type(drop_surface_weights_type), intent(out) :: eff
      end subroutine prepare_surface_weights

      !> [deriv/weights.f90] Copy the cavity's grid-level scalars into a seed state
      !>
      !> @param[in]    self           DROP cavity instance
      !> @param[in]    igrid          Grid point to describe
      !> @param[in]    want_curvature Whether the curvature invariants are needed
      !> @param[inout] state          Seed state
      module subroutine fill_seed_state(self, igrid, want_curvature, state)
         implicit none (type, external)
         class(cavity_type_drop), intent(in) :: self
         integer, intent(in) :: igrid
         logical, intent(in) :: want_curvature
         type(drop_seed_state_type), intent(inout) :: state
      end subroutine fill_seed_state

      !> Contract surface weights to per-grid LSF adjoint weights
      !>
      !> @param[in]  self       DROP cavity instance
      !> @param[in]  acc        Accumulated surface-observable adjoints
      !> @param[out] w_lsf0     Adjoint weights for LSF values (ngrid)
      !> @param[out] w_lsf1     Adjoint weights for LSF gradients (3, ngrid)
      !> @param[out] w_lsf2     Adjoint weights for LSF Hessians (3, 3, ngrid)
      !> @param[out] error      Error object
      module subroutine contract_surface_lsf_weights(self, acc, w_lsf0, w_lsf1, w_lsf2, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(in) :: self
         type(cavity_surface_adjoint_type), intent(in) :: acc
         real(wp), intent(out) :: w_lsf0(:)
         real(wp), intent(out) :: w_lsf1(:, :)
         real(wp), intent(out) :: w_lsf2(:, :, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine contract_surface_lsf_weights

      !> [deriv/potential.f90] Map accumulated surface adjoints into the generic potential
      module subroutine get_surface_potential_drop(self, acc, potential, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(cavity_surface_adjoint_type), intent(in) :: acc
         type(potential_type), intent(inout) :: potential
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_surface_potential_drop

      !> [deriv/nuclear.f90] Contract accumulated surface adjoints into dE/dR_A
      !>
      !> @param[in]    self     DROP cavity instance
      !> @param[in]    acc      Accumulated surface-observable adjoints
      !> @param[inout] gradient Nuclear-gradient accumulator (3, nsph)
      !> @param[out]   error    Error object
      module subroutine get_surface_gradient_drop(self, acc, gradient, error)
         implicit none (type, external)
         class(cavity_type_drop), intent(in) :: self
         type(cavity_surface_adjoint_type), intent(in) :: acc
         real(wp), intent(inout) :: gradient(:, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_surface_gradient_drop

   end interface

contains

   !* ================================================================================= *!
   !*                                    Constructor                                    *!
   !* ================================================================================= *!

   !> Initialize DROP model
   !>
   !> Sets up a DROP cavity instance with optional settings/configuration
   !>
   !> The LSF model is *required*; callers build their LSF concrete (e.g. `svdw%new(...)`)
   !> and pass it as `lsf_model`
   !>
   !> The cavity pushes its derived `screening_threshold` into the LSF so the LSF's internal screening
   !> caches stay consistent with the cavity tolerance
   !>
   !>
   !> @param[inout] self          Cavity instance to initialize
   !> @param[in]    ctx           Shared run context (borrowed; must outlive the cavity)
   !> @param[in]    nleb          Number of Lebedev points per sphere for angular grid (optional)
   !> @param[in]    tolerance     Master numerical tolerance (optional)
   !> @param[in]    proj_maxiter  Maximum number of projection iterations (optional)
   !> @param[in]    proj_level    Projection refinement level (optional)
   !> @param[in]    radius_model  Atomic radius model to use for cavity construction
   !> @param[in]    lsf_model     LSF template (required; cavity stores a copy)
   !> @param[in]    do_fine      Enable all optional properties (optional)
   !> @param[out]   error         Error handling structure (optional)
   !> Initialize DROP cavity
   subroutine new_cavity_drop(self, &
                              ctx, &
                              nleb, &
                              tolerance, proj_maxiter, proj_level, &
                              wleb_prune_level, &
                              do_fine, &
                              radius_model, &
                              lsf_model, &
                              error)
      type(cavity_type_drop), intent(inout) :: self

      !> Shared run context (verbosity/debug/timer); borrowed, must outlive self
      type(moist_context_type), intent(in), target :: ctx

      !> Grid settings
      integer, intent(in), optional :: nleb

      !> Master numerical tolerance
      real(wp), intent(in), optional :: tolerance
      integer, intent(in), optional :: proj_maxiter
      integer, intent(in), optional :: proj_level

      !> Weight switching level (0=off, 1-6=increasing aggressiveness)
      integer, intent(in), optional :: wleb_prune_level

      !> Enable all optional properties
      logical, intent(in), optional :: do_fine

      !> Radius model to use for cavity construction (provided by caller)
      class(radius_type), intent(in) :: radius_model

      !> LSF model template (provided by caller)
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_model

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Borrow the shared run context (owns verbosity/debug/timer)
      self%ctx => ctx

      !> Convenience property shortcuts
      if (present(do_fine)) then
         if (do_fine) self%request = drop_request_fine()
      end if

      !> Parameter setup
      call self%param%new( &
         nleb=nleb, &
         tolerance=tolerance, proj_maxiter=proj_maxiter, proj_level=proj_level, &
         branch_weight_s=0.05_wp, &
         wleb_prune_level=wleb_prune_level, &
         error=error)
      if (allocated(error)) return

      !> Radius model setup
      if (allocated(self%radius_model)) deallocate (self%radius_model)
      allocate (self%radius_model, source=radius_model)

      !> LSF model setup
      if (allocated(self%lsf_model)) deallocate (self%lsf_model)
      allocate (self%lsf_model, source=lsf_model)
      !> Push the cavity-derived screening threshold into the LSF
      self%lsf_model%screening_threshold = self%param%screening_threshold

      !> Set up weight switching function
      call new_sigmoid_bump_swif(self%f_crit, self%param%w_0ls_from, self%param%w_0ls_to, &
                                 p_hi=self%param%w_0ls_p, a_hi=self%param%w_0ls_a, p_lo=self%param%w_0ls_p, a_lo=self%param%w_0ls_a)

      call new_sigmoid_bump_swif(self%f_foc, self%param%w_0tra_from, self%param%w_0tra_to, &
                                 p_hi=self%param%w_0ls_p, a_hi=self%param%w_0ls_a, p_lo=self%param%w_0ls_p, a_lo=self%param%w_0ls_a)

      !> Set up Lebedev weight switching function (optional)
      if (self%param%wleb_prune_level > 0) then
         call new_sigmoid_bump_swif(self%f_wleb, &
                                    self%param%wleb_prune_from, self%param%wleb_prune_to)
      end if

      !> Set up hard-sphere gaussians function (iSwiG)
      call new_iswig(self%iswig, self%param%iswig_xi_born)

      !> Set up branch weight model
      call self%branch_weight%init(self%param%branch_weight_s)

      ! Print parameters and request
      if (self%ctx%verbosity > 1) then
         call self%param%print(unit=self%ctx%unit)
         select type (m => self%lsf_model)
         type is (moist_cavity_drop_lsf_svdw_type)
            call m%param%print(unit=self%ctx%unit)
         end select
         call self%request%print(unit=self%ctx%unit)
      end if

   end subroutine new_cavity_drop

   !> Configure which optional properties to compute and store
   !>
   !> @param[inout] self              Cavity instance
   !> @param[in]    do_fine           Enable all optional properties (optional)
   !> @param[in]    do_curvature      Compute mean and Gaussian curvatures (optional)
   !> @param[in]    do_grid_density   Compute local grid-point density (optional)
   !> @param[in]    do_normal         Store surface normal vectors (optional)
   !> @param[in]    do_r_iI           Store sphere-center to grid-point distances (optional)
   !> @param[in]    do_rho            Store anchor-to-projected-point displacements (optional)
   subroutine set_properties_drop(self, &
                                  do_fine, do_curvature, do_grid_density, &
                                  do_normal, do_r_iI, do_rho)
      class(cavity_type_drop), intent(inout) :: self
      logical, intent(in), optional :: do_fine
      logical, intent(in), optional :: do_curvature
      logical, intent(in), optional :: do_grid_density
      logical, intent(in), optional :: do_normal
      logical, intent(in), optional :: do_r_iI
      logical, intent(in), optional :: do_rho

      !> do_fine sets everything at once
      if (present(do_fine)) then
         if (do_fine) self%request = drop_request_fine()
      end if

      !> Individual flags (applied after do_fine so they can override)
      if (present(do_curvature)) self%request%curvature = do_curvature
      if (present(do_grid_density)) self%request%grid_point_density = do_grid_density
      if (present(do_normal)) self%request%normal = do_normal
      if (present(do_r_iI)) self%request%r_iI = do_r_iI
      if (present(do_rho)) self%request%rho = do_rho

   end subroutine set_properties_drop

   !* ================================================================================= *!
   !*                          Update Cavity (construct cavity)                         *!
   !* ================================================================================= *!

   !> Update cavity for new geometry
   !>
   !> @param[inout] self  Cavity instance
   !> @param[in]    mol   Molecular structure
   subroutine update_cavity_drop(self, mol, error)
      class(cavity_type_drop), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error

      !> Timer stack depth at entry; error paths unwind back to it (below).
      integer :: d0

      !> Set number of spheres
      self%nsph = mol%nat

      call self%radius_model%update(mol, error)
      if (allocated(error)) return
      if (self%ctx%verbosity >= 2) call self%radius_model%print()
      if (allocated(self%radii)) deallocate (self%radii)
      allocate (self%radii(self%nsph))
      self%radii = self%radius_model%f0

      !> Set centers of spheres
      self%mol = mol
      if (allocated(self%sphxyz)) deallocate (self%sphxyz)
      allocate (self%sphxyz(3, self%nsph), source=mol%xyz)

      !> Refresh LSF geometry caches
      call self%lsf_model%update(self%mol, self%radii)

      !* --------------------------------- Setup phase -------------------------------- *!
      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("Setup", category=cat_setup)

      !> Ensure Lebedev cache for current num_leb
      call self%ctx%timer%start("Lebedev cache")
      call self%ensure_lebedev_cache(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Lebedev cache")

      !> Fill intermediate arrays
      call self%ctx%timer%start("Array setup")
      call self%fill_arrays(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Array setup")

      !> Neighbour list setups (before switching so cell grid is available for screening)
      call self%ctx%timer%start("Adj. lists")
      call self%setup_mol_cell_grid(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Adj. lists")

      !> Compute switching function (uses cell grid for O(1) candidate lookup)
      call self%ctx%timer%start("Switching func.")
      call self%compute_switching(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Switching func.")

      !> Pre-filter points below switching cutoff (reduces projection workload)
      call self%ctx%timer%start("Pre-filter")
      call self%filter_arrays('Prefilter', error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Pre-filter")

      call self%ctx%timer%stop("Setup")

      !> Adjust nmax (i.e. the number of anchor points)
      self%nmax = self%ngrid

      !* ------------------------------ Projection phase ------------------------------ *!
      call self%ctx%timer%start("Projector", category=cat_solve)
      call self%project_all_points(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Projector")

      !* ---------------------------- Post processing phase --------------------------- *!
      call self%ctx%timer%start("Post processing", category=cat_setup)

      !> Compute closest-point Jacobian scaling
      call self%ctx%timer%start("CP Jacobian")
      call self%compute_cp_jacobian_scaling(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("CP Jacobian")

      !> Adjust integration weights for scaling and switching
      self%wleb = self%wleb*self%cpjac_scal0*self%w_f0

      !> Switch off small lebedev weights
      if (self%param%wleb_prune_level > 0) then
         block
            integer :: i
            do i = 1, self%ngrid
               self%wleb(i) = self%wleb(i)*self%f_wleb%f0(abs(self%wleb(i)))
            end do
         end block
      end if

      !> Compute branch weights over the final surviving anchor groups
      call self%ctx%timer%start("Branch weights")
      call self%compute_branch_weights(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Branch weights")

      !> Filter out points below cutoff
      call self%ctx%timer%start("Filter")
      call self%filter_arrays('Postfilter', error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Filter")

      !> Setup grid point adjacency list for density computation
      call self%ctx%timer%start("Grid adj. list")
      call self%setup_grid_adj_list(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Grid adj. list")

      !> Nearest neighbour search to find disconnected cavities
      call self%ctx%timer%start("Disconnected cav.")
      call self%find_disconnected_cavities(error=error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Disconnected cav.")

      call self%ctx%timer%stop("Post processing")

      !* ------------------------------ Properties phase ------------------------------ *!
      call self%ctx%timer%start("Properties", category=cat_properties)

      !> Compute grid point densities [optional diagnostic]
      if (self%request%grid_point_density) then
         call self%ctx%timer%start("Grid density")
         call self%compute_grid_point_density(error)
         if (allocated(error)) then
            call self%ctx%timer%unwind(d0)
            return
         end if
         call self%ctx%timer%stop("Grid density")
      end if

      !> Compute curvatures [optional diagnostic]
      if (self%request%curvature) then
         call self%ctx%timer%start("Curvatures")
         call self%compute_curvature(error)
         if (allocated(error)) then
            call self%ctx%timer%unwind(d0)
            return
         end if
         call self%ctx%timer%stop("Curvatures")
      end if

      !> Compute area and volume
      call self%ctx%timer%start("Area & Volume")
      call self%compute_area_volume(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Area & Volume")

      !> Compute Gaussian surface charge widths
      call self%ctx%timer%start("Gaussians")
      call self%compute_gaussians(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Gaussians")

      call self%ctx%timer%stop("Properties")

      !> Run grid diagnostics
      if (self%ctx%verbosity >= 2) then
         call self%ctx%timer%start("Analysis", category=cat_properties)
         call self%analyze_cavity(error)
         if (allocated(error)) then
            call self%ctx%timer%stop("Analysis")
            return
         end if
         call self%ctx%timer%stop("Analysis")
      end if

   end subroutine update_cavity_drop

   !* ================================================================================= *!
   !*                                 First derivatives                                 *!
   !* ================================================================================= *!

   !> Compute and store all requested DROP nuclear derivatives.
   !> @param[inout] self  Cavity instance receiving the derivatives
   !> @param[out]   error Error handling
   subroutine get_gradient_drop(self, error)
      class(cavity_type_drop), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Timer stack depth at entry; error paths unwind back to it (below).
      integer :: d0

      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("Gradients", category=cat_gradient)

      if (self%ctx%verbosity > 1) write (self%ctx%unit, '(a)') "[Info] Computing gradients ..."
      call self%compute_gradient_drop(error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if

      call self%ctx%timer%stop("Gradients")

   end subroutine get_gradient_drop

   !* ================================================================================= *!
   !*                                Second derivatives                                 *!
   !* ================================================================================= *!

   ! Coming soon ;)

   !* ================================================================================= *!
   !*                            Grids and cashing/screening                            *!
   !* ================================================================================= *!

   !> Ensure Lebedev grid cache is initialized and matches the requested size
   ! TODO: A simple wrapper for this into the lebedev grid module would be better
   ! (code deduplication as its also used in iswig and numsa,..)
   subroutine ensure_lebedev_cache(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: oleb, i

      ! Map requested num_leb to Lebedev order index
      call lebedev_order_from_num(self%param%num_leb, oleb, error)
      if (allocated(error)) return

      if (allocated(self%ang_grid) &
          .and. allocated(self%ang_weight) &
          .and. allocated(self%oleb) &
          ) then
         if (.not. allocated(self%nmax)) allocate (self%nmax)
         self%nmax = self%param%num_leb*self%nsph
         return
      end if

      if (allocated(self%ang_grid)) deallocate (self%ang_grid)
      if (allocated(self%ang_weight)) deallocate (self%ang_weight)
      if (allocated(self%oleb)) deallocate (self%oleb)
      if (allocated(self%nmax)) deallocate (self%nmax)

      allocate (self%oleb)
      self%oleb = oleb

      allocate (self%ang_grid(3, self%param%num_leb))
      allocate (self%ang_weight(self%param%num_leb))
      call get_angular_grid(self%oleb, self%ang_grid, self%ang_weight, error)
      if (allocated(error)) return

      !> Check for negative weights (?!)
      if (any(self%ang_weight < 0.0_wp)) then
         call fatal_error(error, "Grid contains negativ weights that do not work with DROP.")
         return
      end if

      allocate (self%nmax)
      self%nmax = self%param%num_leb*self%nsph

   end subroutine ensure_lebedev_cache

   !* ================================================================================= *!
   !*                                       Debug                                       *!
   !* ================================================================================= *!

   !> Write grid to csv file (for debugging)
   subroutine write_cavity_csv_debug(self, filename, error)
      class(cavity_type_drop), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, i
      real(wp) :: val_wleb, val_r_iI0, val_rho
      real(wp) :: val_anch_x, val_anch_y, val_anch_z
      real(wp) :: val_n_x, val_n_y, val_n_z
      real(wp) :: val_rho_grid, val_rho_grid_anchor, val_KM, val_KG, val_cpjac
      real(wp) :: val_sigma_max, val_sigma_min, val_sigma_chi
      real(wp) :: val_risk_score, val_gamma_tilde, val_det_B, val_kappa_B
      logical :: val_converged

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, 'write_csv_debug: cavity grid not allocated')
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, 'write_csv_debug: no grid points to write')
         return
      end if

      open (file=filename, newunit=unit, status='replace', action='write', iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, 'Could not open CSV file for writing: '//trim(filename))
         return
      end if

      write (unit, '(a)') 'ngrid,numbering,x,y,z,owner,area,switch_f,w_leb,rho_scal0,rad,r_iI0,rho,'// &
         'anch_x,anch_y,anch_z,n_x,n_y,n_z,rho_grid,rho_grid_anchor,KM,KG,cpjac_scal,'// &
         ',converged'

      do i = 1, self%ngrid

         val_rho_grid = 0.0_wp
         val_rho_grid_anchor = 0.0_wp
         val_KM = 0.0_wp
         val_KG = 0.0_wp
         val_cpjac = 1.0_wp

         ! Safely extract values
         if (allocated(self%rho_grid)) val_rho_grid = self%rho_grid(i)
         if (allocated(self%rho_grid_anchor)) val_rho_grid_anchor = self%rho_grid_anchor(i)
         if (allocated(self%k1)) val_KM = self%k1(i)
         if (allocated(self%k2)) val_KG = self%k2(i)
         if (allocated(self%cpjac_scal0)) val_cpjac = self%cpjac_scal0(i)
         if (allocated(self%w_f0)) val_kappa_B = self%w_f0(i)

         write (unit, '(i0,",",i0,30(",",g0))') &
            i, &
            self%numbering(i), &
            self%xyz(1, i), self%xyz(2, i), self%xyz(3, i), &
            self%owner(i), &
            self%a(i), &
            self%f(i), &
            self%wleb(i), 0.0_wp, &
            self%radii(self%owner(i)), &
            self%r_iI0(i), &
            self%rho(i), &
            self%anchorxyz(1, i), self%anchorxyz(2, i), self%anchorxyz(3, i), &
            self%normal0(1, i), self%normal0(2, i), self%normal0(3, i), &
            val_rho_grid, val_rho_grid_anchor, &
            val_KM, val_KG, val_cpjac, &
            self%converged(i)
      end do
      close (unit)

      write (self%ctx%unit, '(a,1x,a)') '[Info] Wrote cavity grid to', trim(filename)

   end subroutine write_cavity_csv_debug

   !* ================================================================================= *!
   !*                                     Finalizer                                     *!
   !* ================================================================================= *!

   !> Finalizer for cavity_type_drop to properly deallocate all allocatable components
   !> This ensures proper cleanup when the cavity is deleted through the C API
   subroutine finalize_cavity_drop(self)
      type(cavity_type_drop), intent(inout) :: self

      !> The profiling timer lives on the borrowed run context (self%ctx) and is
      !> owned by the top-level caller, so it is not torn down here. Callers that
      !> want the accumulated tree call self%ctx%timer%write(...) before teardown.

      ! Deallocate grid point data arrays
      if (allocated(self%xi0)) deallocate (self%xi0)
      if (allocated(self%f)) deallocate (self%f)
      if (allocated(self%wleb)) deallocate (self%wleb)
      if (allocated(self%anchor_wleb0)) deallocate (self%anchor_wleb0)
      if (allocated(self%rho)) deallocate (self%rho)
      if (allocated(self%r_iI0)) deallocate (self%r_iI0)
      if (allocated(self%normal0)) deallocate (self%normal0)
      if (allocated(self%anchorxyz)) deallocate (self%anchorxyz)
      if (allocated(self%numbering)) deallocate (self%numbering)
      if (allocated(self%v)) deallocate (self%v)
      if (allocated(self%vsph)) deallocate (self%vsph)

      ! Deallocate grid settings/caching
      if (allocated(self%nmax)) deallocate (self%nmax)
      if (allocated(self%oleb)) deallocate (self%oleb)
      if (allocated(self%ang_grid)) deallocate (self%ang_grid)
      if (allocated(self%ang_weight)) deallocate (self%ang_weight)

      ! Deallocate gradient arrays
      if (allocated(self%lambda0)) deallocate (self%lambda0)
      if (allocated(self%branch)) deallocate (self%branch)
      if (allocated(self%anchor_id)) deallocate (self%anchor_id)
      if (allocated(self%branch_count)) deallocate (self%branch_count)
      if (allocated(self%wbranch)) deallocate (self%wbranch)
      if (allocated(self%phi0)) deallocate (self%phi0)
      if (allocated(self%xyz1_rA)) deallocate (self%xyz1_rA)
      if (allocated(self%xyz2_rArB)) deallocate (self%xyz2_rArB)
      if (allocated(self%rho1_rA)) deallocate (self%rho1_rA)
      if (allocated(self%lambda1_rA)) deallocate (self%lambda1_rA)
      if (allocated(self%lambda2_rArB)) deallocate (self%lambda2_rArB)
      if (allocated(self%r_iI1_rA)) deallocate (self%r_iI1_rA)
      if (allocated(self%a_i1_rA)) deallocate (self%a_i1_rA)
      if (allocated(self%v1_rA)) deallocate (self%v1_rA)
      if (allocated(self%xi1_rA)) deallocate (self%xi1_rA)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      if (allocated(self%cpjac_scal0)) deallocate (self%cpjac_scal0)
      if (allocated(self%A_tot1_rA)) deallocate (self%A_tot1_rA)
      if (allocated(self%asph1_rA)) deallocate (self%asph1_rA)
      if (allocated(self%V_tot1_rA)) deallocate (self%V_tot1_rA)
      if (allocated(self%vsph1_rA)) deallocate (self%vsph1_rA)

      ! Deallocate neighbour list
      call self%mol_cell_grid%destroy()
      call self%grid_adj_list%destroy()

      ! Deallocate structure_type allocatable components
      if (allocated(self%mol%id)) deallocate (self%mol%id)
      if (allocated(self%mol%num)) deallocate (self%mol%num)
      if (allocated(self%mol%sym)) deallocate (self%mol%sym)
      if (allocated(self%mol%xyz)) deallocate (self%mol%xyz)
      if (allocated(self%mol%lattice)) deallocate (self%mol%lattice)
      if (allocated(self%mol%periodic)) deallocate (self%mol%periodic)
      if (allocated(self%mol%bond)) deallocate (self%mol%bond)
      if (allocated(self%mol%comment)) deallocate (self%mol%comment)
      if (allocated(self%mol%sdf)) deallocate (self%mol%sdf)
      if (allocated(self%mol%pdb)) deallocate (self%mol%pdb)

      ! Deallocate inherited allocatable components from base cavity_type
      if (allocated(self%asph)) deallocate (self%asph)
      if (allocated(self%total_area)) deallocate (self%total_area)
      if (allocated(self%total_volume)) deallocate (self%total_volume)
      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%sphxyz)) deallocate (self%sphxyz)
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%a)) deallocate (self%a)
      if (allocated(self%owner)) deallocate (self%owner)
      if (allocated(self%converged)) deallocate (self%converged)

   end subroutine finalize_cavity_drop

end module moist_cavity_drop
