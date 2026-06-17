!> Main DROP (Discretization via Reference-Onto surface Projection) implementation
module moist_cavity_drop
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_structure, only: structure_type
   use mctc_io, only: new
   use mctc_env, only: error_type, fatal_error, wp
   use iso_fortran_env, only: error_unit, output_unit
   use moist_math_lapack, only: getrf, getrs
   use moist_math_lapack_gesv, only: dgesv
   use moist_math_linalg, only: mat3x3_inv, setup_tangent_frame
   use moist_math_boys, only: dboysfun1
   use moist_math_grid_lebedev, only: get_angular_grid, grid_size, lebedev_order_from_num
   use moist_type, only: cavity_type
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
   use moist_utils_prettyprint, only: prettyprinter
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter

   use moist_cavity_drop_marchingcubes, only: integrate_surface_marching_cubes

   use moist_math_smoothing_kernels, only: wendland_kernel_type

   use moist_utils_timer, only: timer_type
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

      !> Level-set function model. Constructed once at cavity setup; every
      !> per-thread LSF (in projector, projection, gradient, properties,
      !> marching cubes) is sourced-allocated from this template.
      class(moist_cavity_drop_lsf_type), allocatable :: lsf_model

      !> Property request flags controlling which quantities are computed
      type(drop_property_request) :: request

      !> Timer for profiling cavity update steps
      type(timer_type) :: timer

      !> ======== Atomic sphere data ========

      !> Molecular structure
      type(structure_type) :: mol

      !> Number of atomic spheres
      integer :: nsph

      !> Per-cell candidate atom lists for point-to-atom screening
      type(moist_cell_grid_type) :: mol_cell_grid

      !> Grid-point neighbour list for density computation (CSR format)
      type(adjacency_list_type) :: grid_adj_list

      !> Unique numbering
      integer, allocatable :: numbering(:)

      !> ======== Integration grid ========

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

      !> ======== Geometric quantities ========

      !> Initial (unprojected) grid positions
      real(wp), allocatable :: anchorxyz(:, :)

      !> Per-gridpoint derivatives of gridpoint w.r.t. nuclear coordinates
      !> Shape: (3, 3, nsph, ngrid) = (j, alpha, A, igrid)
      real(wp), allocatable :: xyz1_rA(:, :, :, :)
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

      !> Surface normal vectors at grid points
      real(wp), allocatable :: normal0(:, :)
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

      !> ======== Areas ========

      !> Per-gridpoint derivatives of the area (3, nsph, ngrid)
      real(wp), allocatable :: a_i1_rA(:, :, :)
      !> Accumulated gradient of per-atom areas w.r.t. nuclear coordinates (3, nsph, nsph_owner)
      real(wp), allocatable :: asph1_rA(:, :, :)
      !> Accumulated gradient of total area w.r.t. nuclear coordinates (3, nsph)
      real(wp), allocatable :: A_tot1_rA(:, :)

      !> ======== Volumes ========

      !> Grid point volumes
      real(wp), allocatable :: v(:)
      !> Per-gridpoint derivatives of the volume (3, nsph, ngrid)
      real(wp), allocatable :: v_i1_rA(:, :, :)
      !> Atomic volumes
      real(wp), allocatable :: vsph(:)
      !> Accumulated gradient of per-atom volumes w.r.t. nuclear coordinates (3, nsph, nsph)
      real(wp), allocatable :: vsph1_rA(:, :, :)
      !> Accumulated gradient of total volume w.r.t. nuclear coordinates (3, nsph)
      real(wp), allocatable :: V_tot1_rA(:, :)

      !> ======== Grid settings and caching ========

      !> Total number of initial grid points (before filtering)
      integer, allocatable :: nmax
      !> Lebedev grid order index
      integer, allocatable :: oleb
      !> Cached Lebedev angular grid (3, num_leb)
      real(wp), allocatable :: ang_grid(:, :)
      !> Cached Lebedev weights (num_leb)
      real(wp), allocatable :: ang_weight(:)

      !> ======== Switching functions ========

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

      !> Combined switching function value (ngrid)
      real(wp), allocatable :: f(:)
      !> Per-gridpoint derivatives of switching function f w.r.t. nuclear coordinates (3, nsph, ngrid)
      real(wp), allocatable :: f1_rA(:, :, :)

      !> iSwiG switching values per point
      real(wp), allocatable :: iswig_f0(:)

      !> ======== Projection ========

      !> Branch weight model for handling degenerate branches
      type(branch_weight_type) :: branch_weight

      !> Branch index per grid point (1 = default/unbranched)
      integer, allocatable :: branch(:)
      !> Anchor group id for each grid point
      integer, allocatable :: anchor_id(:)
      !> Number of branches in the anchor group of each grid point
      integer, allocatable :: branch_count(:)
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

      !> Verbosity level (0=quiet, 1=normal, 2=verbose)
      integer :: verbosity = 2
      logical :: debug = .false.

      !> ======== CPCM ========

      !> Gaussian width at each anchor grid point (ngrid)
      real(wp), allocatable :: anchor_xi0(:)
      !> Gaussian width at each grid point (ngrid)
      real(wp), allocatable :: xi0(:)
      !> Per-gridpoint derivatives of xi (3, nsph, ngrid)
      real(wp), allocatable :: xi1_rA(:, :, :)

      !> ======== CPCM debug/temporary storage ========
      ! TODO: These routines are/were used for development/testing and will be removed at some point

      !> CPCM surface charges (ngrid) [DEBUG/TEMPORARY]
      real(wp), allocatable :: cpcm_q(:)
      !> CPCM electrostatic potential at grid points (ngrid) [DEBUG/TEMPORARY]
      real(wp), allocatable :: cpcm_pot(:)
      !> Fixed source charges used by CPCM test model (nsph) [DEBUG/TEMPORARY]
      real(wp), allocatable :: cpcm_source_charges(:)
      !> CPCM solvation energy [DEBUG/TEMPORARY]
      real(wp) :: cpcm_energy
      !> CPCM energy gradient (3, nsph) [DEBUG/TEMPORARY]
      real(wp), allocatable :: cpcm_gradient(:, :)

   contains
      !> Configure which optional properties to compute
      procedure :: properties => set_properties_drop
      !> Update cavity for new geometry
      procedure :: update => update_cavity_drop
      !> Compute area gradient w.r.t. nuclear coordinates
      procedure :: get_gradient => get_gradient_drop

      !> Compute CPCM properties
      procedure :: compute_gaussians

      !> Assemble PCM interaction matrix (DROP Gaussian-based A-matrix)
      procedure :: get_amat => get_amat_drop
      !> Compute CPCM derivatives
      procedure :: Amat012_rA => assemble_Amat012_rA
      !> Contract CPCM A-matrix first derivatives with two grid vectors
      procedure :: contract_amat1_q1q2_rA
      !> Contract CPCM A-matrix derivatives to per-grid surface weights
      procedure :: contract_amat1_q1q2_surface_weights
      !> Contract surface-coordinate weights to LSF adjoint weights
      procedure :: contract_surface_lsf_weights
      !> Contract combined nuclear and electronic terms for CPCM gradients
      procedure :: contract_nuc_elec_qefield_rA

      !> Compute all needed cavity gradients
      procedure :: compute_gradient_drop

      !> Compute CPCM solvation energy (debug/temporary)
      procedure :: compute_cpcm_energy
      !> Compute CPCM energy gradient (debug/temporary)
      procedure :: compute_cpcm_energy_gradient

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
         integer, intent(in) :: current_capacity
         integer, intent(in) :: required_capacity
      end function projection_grow_capacity

      !> [projection.f90] Ensure projection arrays are allocated for `new_capacity`
      module subroutine ensure_projection_capacity(self, new_capacity)
         class(cavity_type_drop), intent(inout) :: self
         integer, intent(in) :: new_capacity
      end subroutine ensure_projection_capacity

      !> [projection.f90] Project all current grid points onto the LSF surface
      module subroutine project_all_points(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine project_all_points

      !> [setup.f90] Evaluate switching weights for all grid points
      module subroutine compute_switching_function(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_switching_function

      !> [properties.f90] Compute local grid-point density values
      module subroutine compute_grid_point_density(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_grid_point_density

      !> [properties.f90] Check cavity diagnostics
      module subroutine analyze_cavity(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine analyze_cavity

      !> [setup.f90] Fill initial DROP arrays from per-atom Lebedev grids
      module subroutine fill_arrays(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine fill_arrays

      !> [setup.f90] Build per-cell atom screening grid
      module subroutine setup_mol_cell_grid(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine setup_mol_cell_grid

      !> [setup.f90] Build grid-point neighbour list for density computation
      module subroutine setup_grid_adj_list(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine setup_grid_adj_list

      !> [filter.f90] Filter and compact points after projection
      module subroutine filter_arrays(self, name, error)
         class(cavity_type_drop), intent(inout) :: self
         character(len=*), intent(in) :: name
         type(error_type), allocatable, intent(out) :: error
      end subroutine filter_arrays

      !> [filter.f90] Compute branch weights on the final surviving set
      module subroutine compute_branch_weights(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_branch_weights

      !> [properties.f90] Compute total/atomic area and volume contributions
      module subroutine compute_area_volume(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_area_volume

      !> [properties.f90] Compute mean and Gaussian curvature on the surface grid
      module subroutine compute_curvature(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_curvature

      !> [projection.f90] Compute closest-point Jacobian scaling factors.
      module subroutine compute_cp_jacobian_scaling(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_cp_jacobian_scaling

      !> [projection.f90] Compute CPCM Gaussian width parameters for grid points
      module subroutine compute_gaussians(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_gaussians

      !> [gradient.f90] Compute first nuclear derivatives for DROP quantities
      module subroutine compute_gradient_drop(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_gradient_drop

      !> [cpcm.f90] Wrapper around assemble_Amat012_rA for the cavity get_amat interface
      module subroutine get_amat_drop(self, amat, error)
         class(cavity_type_drop), intent(in) :: self
         real(wp), intent(out) :: amat(:, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_amat_drop

      !* ============================================================================== *!
      !*                         CPCM routines for QM interface                         *!
      !* ============================================================================== *!

      module subroutine assemble_Amat012_rA(self, Amat0, Amat1_rA, error)
         class(cavity_type_drop), intent(in) :: self
         real(wp), allocatable, intent(out) :: Amat0(:, :)
         real(wp), allocatable, optional, intent(out) :: Amat1_rA(:, :, :, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine assemble_Amat012_rA

      !> Contract first derivatives of CPCM A with two grid vectors
      !>
      !> Computes:
      !> - `grad_rA = \sum_{ij} q1_i (\partial A_{ij}/\partial R_A) q2_j`
      !>
      !> @param[in]  self        DROP cavity instance
      !> @param[in]  q1          Left contraction vector (ngrid)
      !> @param[in]  q2          Right contraction vector (ngrid)
      !> @param[out] grad_rA     Contracted gradient contribution (3, nsph)
      !> @param[out] error       Error object (optional)
      module subroutine contract_amat1_q1q2_rA(self, q1, q2, grad_rA, error)
         class(cavity_type_drop), intent(in) :: self
         real(wp), intent(in) :: q1(:)
         real(wp), intent(in) :: q2(:)
         real(wp), intent(out) :: grad_rA(3, self%nsph)
         type(error_type), allocatable, intent(out) :: error
      end subroutine contract_amat1_q1q2_rA

      !> Contract first derivatives of CPCM A to per-grid surface weights
      !>
      !> @param[in]  self    DROP cavity instance
      !> @param[in]  q1      Left contraction vector (ngrid)
      !> @param[in]  q2      Right contraction vector (ngrid)
      !> @param[out] w_xi    Weights for Gaussian widths `xi_i` (ngrid)
      !> @param[out] w_f     Weights for switch factors `f_i` (ngrid)
      !> @param[out] w_xyz   Weights for surface coordinates `r_i` (3, ngrid)
      !> @param[out] error   Error object
      module subroutine contract_amat1_q1q2_surface_weights(self, q1, q2, w_xi, w_f, w_xyz, error)
         class(cavity_type_drop), intent(in) :: self
         real(wp), intent(in) :: q1(:)
         real(wp), intent(in) :: q2(:)
         real(wp), intent(out) :: w_xi(:)
         real(wp), intent(out) :: w_f(:)
         real(wp), intent(out) :: w_xyz(:, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine contract_amat1_q1q2_surface_weights

      !> Contract surface weights to per-grid LSF adjoint weights
      !>
      !> @param[in]  self       DROP cavity instance
      !> @param[in]  w_xi       Surface weights for Gaussian widths
      !> @param[in]  w_f        Surface weights for anchor switch factors
      !> @param[in]  w_xyz      Surface coordinate weights (3, ngrid)
      !> @param[out] w_lsf0     Adjoint weights for LSF values (ngrid)
      !> @param[out] w_lsf1     Adjoint weights for LSF gradients (3, ngrid)
      !> @param[out] w_lsf2     Adjoint weights for LSF Hessians (3, 3, ngrid)
      !> @param[out] error      Error object
      module subroutine contract_surface_lsf_weights(self, w_xi, w_f, w_xyz, &
            & w_lsf0, w_lsf1, w_lsf2, error)
         class(cavity_type_drop), intent(in) :: self
         real(wp), intent(in) :: w_xi(:)
         real(wp), intent(in) :: w_f(:)
         real(wp), intent(in) :: w_xyz(:, :)
         real(wp), intent(out) :: w_lsf0(:)
         real(wp), intent(out) :: w_lsf1(:, :)
         real(wp), intent(out) :: w_lsf2(:, :, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine contract_surface_lsf_weights

      !> Contract combined nuclear and electronic contributions for CPCM gradients
      !>
      !> @param[in]  self       DROP cavity instance
      !> @param[in]  surface_q  Surface charges `q_i` (ngrid)
      !> @param[in]  qefield    Electronic contribution `Q_i E_elec(i)` (3, ngrid)
      !> @param[in]  za         Nuclear charges `Z_K` (nsph)
      !> @param[out] grad_rA    Contracted gradient contribution (3, nsph)
      !> @param[out] error      Error object (optional)
      module subroutine contract_nuc_elec_qefield_rA(self, surface_q, qefield, za, grad_rA, error)
         class(cavity_type_drop), intent(in) :: self
         real(wp), intent(in) :: surface_q(:)
         real(wp), intent(in) :: qefield(:, :)
         real(wp), intent(in) :: za(:)
         real(wp), intent(out) :: grad_rA(3, self%nsph)
         type(error_type), allocatable, intent(out) :: error
      end subroutine contract_nuc_elec_qefield_rA

      !* ============================================================================== *!
      !*            Debug/testing routines (use random fixed partial charges)            *!
      !* ============================================================================== *!

      module subroutine compute_cpcm_energy(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_cpcm_energy

      module subroutine compute_cpcm_energy_gradient(self, error)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_cpcm_energy_gradient
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
   !> @param[in]    verbose       Verbosity level for output (optional)
   !> @param[in]    debug         Enable debug output (optional)
   !> @param[in]    nleb          Number of Lebedev points per sphere for angular grid (optional)
   !> @param[in]    tolerance     Master numerical tolerance (optional)
   !> @param[in]    proj_maxiter  Maximum number of projection iterations (optional)
   !> @param[in]    proj_level    Projection refinement level (optional)
   !> @param[in]    radius_model  Atomic radius model to use for cavity construction
   !> @param[in]    lsf_model     LSF template (required; cavity stores a copy)
   !> @param[in]    do_cpcm      Enable CPCM solvation energy computation (optional)
   !> @param[in]    do_fine      Enable all optional properties (optional)
   !> @param[out]   error         Error handling structure (optional)
   !> Initialize DROP cavity
   subroutine new_cavity_drop(self, &
                              verbose, debug, &
                              nleb, &
                              tolerance, proj_maxiter, proj_level, &
                              wleb_prune_level, &
                              do_cpcm, do_fine, &
                              radius_model, &
                              lsf_model, &
                              error)
      type(cavity_type_drop), intent(inout) :: self

      !> Printing
      integer, intent(in), optional :: verbose
      logical, intent(in), optional :: debug

      !> Grid settings
      integer, intent(in), optional :: nleb

      !> Master numerical tolerance
      real(wp), intent(in), optional :: tolerance
      integer, intent(in), optional :: proj_maxiter
      integer, intent(in), optional :: proj_level

      !> Weight switching level (0=off, 1-6=increasing aggressiveness)
      integer, intent(in), optional :: wleb_prune_level

      !> Enable CPCM solvation energy computation
      logical, intent(in), optional :: do_cpcm
      !> Enable all optional properties
      logical, intent(in), optional :: do_fine

      !> Radius model to use for cavity construction (provided by caller)
      class(radius_type), intent(in) :: radius_model

      !> LSF model template (provided by caller)
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_model

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Debug & Verbosity
      if (present(debug)) self%debug = debug
      if (present(verbose)) self%verbosity = verbose

      !> Convenience property shortcuts
      if (present(do_fine)) then
         if (do_fine) self%request = drop_request_fine()
      end if
      if (present(do_cpcm)) then
         if (do_cpcm) self%request%cpcm = .true.
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
      if (self%verbosity > 1) then
         call self%param%print()
         select type (m => self%lsf_model)
         type is (moist_cavity_drop_lsf_svdw_type)
            call m%param%print()
         end select
         call self%request%print()
      end if

      !> Initialize timer with hierarchical structure
      if (self%verbosity > 1) then
         call self%timer%new(32, verbose=.true.)
      else
         call self%timer%new(32, verbose=.false.)
      end if

      !> Register timers
      call self%timer%register(1, 'Setup')
      call self%timer%register(2, 'Lebedev cache', parent=1)
      call self%timer%register(3, 'Array setup', parent=1)
      call self%timer%register(4, 'Adj. lists', parent=1)
      call self%timer%register(5, 'Switching func.', parent=1)
      call self%timer%register(6, 'Pre-filter', parent=1)

      call self%timer%register(7, 'Projector')

      call self%timer%register(8, 'Post processing')
      call self%timer%register(9, 'Filter', parent=8)
      call self%timer%register(10, 'Grid adj. list', parent=8)
      call self%timer%register(11, 'CP Jacobian', parent=8)
      call self%timer%register(32, 'Branch weights', parent=8)
      call self%timer%register(12, 'Disconnected cav.', parent=8)

      call self%timer%register(13, 'Properties')
      call self%timer%register(14, 'Grid density', parent=13)
      call self%timer%register(15, 'Curvatures', parent=13)
      call self%timer%register(16, 'Area & Volume', parent=13)
      call self%timer%register(17, 'Gaussians', parent=13)
      call self%timer%register(18, 'CPCM energy', parent=13)

      call self%timer%register(19, 'Gradients')
      call self%timer%register(20, 'Primitives', parent=19)
      call self%timer%register(21, 'Positions', parent=19)
      call self%timer%register(22, 'Displacement', parent=19)
      call self%timer%register(23, 'Distances', parent=19)
      call self%timer%register(25, 'CP Jacobian', parent=19)
      call self%timer%register(26, 'Gaussian widths', parent=19)
      call self%timer%register(27, 'Switching func.', parent=19)
      call self%timer%register(28, 'Area', parent=19)
      call self%timer%register(29, 'Volume', parent=19)
      call self%timer%register(30, 'Surface normal', parent=19)
      call self%timer%register(31, 'Branch weights', parent=19)
      call self%timer%register(33, 'CPCM', parent=19)

   end subroutine new_cavity_drop

   !> Configure which optional properties to compute and store
   !>
   !> @param[inout] self              Cavity instance
   !> @param[in]    do_fine           Enable all optional properties (optional)
   !> @param[in]    do_cpcm           Compute CPCM solvation energy (optional)
   !> @param[in]    do_curvature      Compute mean and Gaussian curvatures (optional)
   !> @param[in]    do_grid_density   Compute local grid-point density (optional)
   !> @param[in]    do_normal         Store surface normal vectors (optional)
   !> @param[in]    do_r_iI           Store sphere-center to grid-point distances (optional)
   !> @param[in]    do_rho            Store anchor-to-projected-point displacements (optional)
   !> @param[in]    do_mc             Compute marching-cubes reference area/volume (optional)
   subroutine set_properties_drop(self, &
                                  do_fine, do_cpcm, do_curvature, do_grid_density, &
                                  do_normal, do_r_iI, do_rho, do_mc)
      class(cavity_type_drop), intent(inout) :: self
      logical, intent(in), optional :: do_fine
      logical, intent(in), optional :: do_cpcm
      logical, intent(in), optional :: do_curvature
      logical, intent(in), optional :: do_grid_density
      logical, intent(in), optional :: do_normal
      logical, intent(in), optional :: do_r_iI
      logical, intent(in), optional :: do_rho
      logical, intent(in), optional :: do_mc

      !> do_fine sets everything at once
      if (present(do_fine)) then
         if (do_fine) self%request = drop_request_fine()
      end if

      !> Individual flags (applied after do_fine so they can override)
      if (present(do_cpcm)) self%request%cpcm = do_cpcm
      if (present(do_curvature)) self%request%curvature = do_curvature
      if (present(do_grid_density)) self%request%grid_point_density = do_grid_density
      if (present(do_normal)) self%request%normal = do_normal
      if (present(do_r_iI)) self%request%r_iI = do_r_iI
      if (present(do_rho)) self%request%rho = do_rho
      if (present(do_mc)) self%request%mc = do_mc

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

      real(wp) :: mc_area, mc_volume

      !> Set number of spheres
      self%nsph = mol%nat

      call self%radius_model%update(mol, error)
      if (allocated(error)) return
      if (self%verbosity >= 2) call self%radius_model%print()
      if (allocated(self%radii)) deallocate (self%radii)
      allocate (self%radii(self%nsph))
      self%radii = self%radius_model%f0

      !> Set centers of spheres
      self%mol = mol

      !> Refresh LSF geometry caches
      call self%lsf_model%update(self%mol, self%radii)

      !* --------------------------------- Setup phase -------------------------------- *!
      call self%timer%measure(1)

      !> Ensure Lebedev cache for current num_leb
      call self%timer%measure(2)
      call self%ensure_lebedev_cache(error)
      if (allocated(error)) return
      call self%timer%measure(2)

      !> Fill intermediate arrays
      call self%timer%measure(3)
      call self%fill_arrays(error)
      if (allocated(error)) return
      call self%timer%measure(3)

      !> Neighbour list setups (before switching so cell grid is available for screening)
      call self%timer%measure(4)
      call self%setup_mol_cell_grid(error)
      if (allocated(error)) return
      call self%timer%measure(4)

      !> Compute switching function (uses cell grid for O(1) candidate lookup)
      call self%timer%measure(5)
      call self%compute_switching(error)
      if (allocated(error)) return
      call self%timer%measure(5)

      !> Pre-filter points below switching cutoff (reduces projection workload)
      call self%timer%measure(6)
      call self%filter_arrays('Prefilter', error)
      if (allocated(error)) return
      call self%timer%measure(6)

      call self%timer%measure(1)

      !> Adjust nmax (i.e. the number of anchor points)
      self%nmax = self%ngrid

      !* ------------------------------ Projection phase ------------------------------ *!
      call self%timer%measure(7)
      call self%project_all_points(error)
      if (allocated(error)) return
      call self%timer%measure(7)

      !* ---------------------------- Post processing phase --------------------------- *!
      call self%timer%measure(8)

      !> Compute closest-point Jacobian scaling
      call self%timer%measure(11)
      call self%compute_cp_jacobian_scaling(error)
      if (allocated(error)) return
      call self%timer%measure(11)

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
      call self%timer%measure(32)
      call self%compute_branch_weights(error)
      if (allocated(error)) return
      call self%timer%measure(32)

      !> Filter out points below cutoff
      call self%timer%measure(9)
      call self%filter_arrays('Postfilter', error)
      if (allocated(error)) return
      call self%timer%measure(9)

      !> Setup grid point adjacency list for density computation
      call self%timer%measure(10)
      call self%setup_grid_adj_list(error)
      if (allocated(error)) return
      call self%timer%measure(10)

      !> Nearest neighbour search to find disconnected cavities
      call self%timer%measure(12)
      call self%find_disconnected_cavities(error=error)
      if (allocated(error)) return
      call self%timer%measure(12)

      call self%timer%measure(8)

      !* ------------------------------ Properties phase ------------------------------ *!
      call self%timer%measure(13)

      !> Compute grid point densities [optional diagnostic]
      if (self%request%grid_point_density) then
         call self%timer%measure(14)
         call self%compute_grid_point_density(error)
         if (allocated(error)) return
         call self%timer%measure(14)
      end if

      !> Compute curvatures [optional diagnostic]
      if (self%request%curvature) then
         call self%timer%measure(15)
         call self%compute_curvature(error)
         if (allocated(error)) return
         call self%timer%measure(15)
      end if

      !> Compute area and volume
      call self%timer%measure(16)
      call self%compute_area_volume(error)
      if (allocated(error)) return
      call self%timer%measure(16)

      !> Compute Gaussian surface charge widths
      call self%timer%measure(17)
      call self%compute_gaussians(error)
      if (allocated(error)) return
      call self%timer%measure(17)

      !> Compute CPCM solvation energy (debug/testing)
      if (self%request%cpcm) then
         call self%timer%measure(18)
         if (self%verbosity > 1) write (output_unit, '(a)') "[Info] Computing CPCM energy ..."
         call self%compute_cpcm_energy(error)
         if (allocated(error)) return
         call self%timer%measure(18)
      end if

      call self%timer%measure(13)

      !> Run grid diagnostics
      if (self%verbosity >= 2) then
         call self%analyze_cavity(error)
         if (allocated(error)) return
      end if

      !> Compute the area and volume using marching cubes
      ! TODO: Marching cubes should be a standalone "cavity model" which uses a given LSF
      if (self%request%mc) then
         block
            class(moist_cavity_drop_lsf_type), allocatable :: lsf
            type(prettylistprinter) :: plp
            real(wp) :: drop_area, drop_vol, da, dv

            allocate (lsf, source=self%lsf_model)
            call lsf%update(mol, self%radius_model%f0)
            call integrate_surface_marching_cubes(lsf, mol%xyz, &
                                                  mc_area, mc_volume, &
                                                  verbosity=self%verbosity, debug=self%debug, &
                                                  target_spacing=0.1_wp)

            if (self%verbosity >= 2) then

               write (output_unit, '(a)') 'Maching cubes:'
               write (output_unit, '(2x,a,t23,a,f20.10,1x,a)') 'Total MC area', ' ... ', &
                  mc_area, 'bohr^2'
               write (output_unit, '(2x,a,t23,a,f20.10,1x,a)') 'Total MC volume', ' ... ', &
                  mc_volume, 'bohr^3'
               write (output_unit, '(2x,a,t23,a,f20.10,1x,a)') 'Spherical MC radius', ' ... ', &
                  (3.0_wp/(4.0_wp*pi)*mc_volume)**(1.0_wp/3.0_wp), 'bohr'

               drop_area = self%total_area
               drop_vol = self%total_volume
               da = mc_area - drop_area
               dv = mc_volume - drop_vol

               plp = new_prettylistprinter( &
                     [16, 16, 16, 16, 16], &
                     [character(16) :: '', 'DROP', 'MC', 'Abs. diff.', 'Rel. diff. (%)'], &
                     unit=output_unit)
               call plp%blank()
               call plp%print_header()
               call plp%separator()

               call plp%add('Area (bohr^2)')
               call plp%add(drop_area)
               call plp%add(mc_area)
               call plp%add(da)
               if (abs(mc_area) > 0.0_wp) then
                  call plp%add(da/mc_area*100.0_wp)
               else
                  call plp%add('N/A')
               end if
               call plp%end_row()

               call plp%add('Volume (bohr^3)')
               call plp%add(drop_vol)
               call plp%add(mc_volume)
               call plp%add(dv)
               if (abs(mc_volume) > 0.0_wp) then
                  call plp%add(dv/mc_volume*100.0_wp)
               else
                  call plp%add('N/A')
               end if
               call plp%end_row()

               call plp%separator()
               call plp%blank()
            end if
         end block
      end if

   end subroutine update_cavity_drop

   !* ================================================================================= *!
   !*                                 First derivatives                                 *!
   !* ================================================================================= *!

   subroutine get_gradient_drop(self)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable :: gradient_error

      call self%timer%measure(19)

      if (self%verbosity > 1) write (output_unit, '(a)') "[Info] Computing gradients ..."
      call self%compute_gradient_drop(gradient_error)
      if (allocated(gradient_error)) return

      !> Compute Amat gradients and CPCM energy
      call self%timer%measure(33)
      if (self%request%cpcm) then
         if (self%verbosity > 1) write (output_unit, '(a)') "[Info] CPCM gradients ..."
         call self%compute_cpcm_energy_gradient(gradient_error)
         if (allocated(gradient_error)) return
      end if
      call self%timer%measure(33)

      call self%timer%measure(19)

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

      write (output_unit, '(a,1x,a)') '[Info] Wrote cavity grid to', trim(filename)

   end subroutine write_cavity_csv_debug

   !* ================================================================================= *!
   !*                                     Finalizer                                     *!
   !* ================================================================================= *!

   !> Finalizer for cavity_type_drop to properly deallocate all allocatable components
   !> This ensures proper cleanup when the cavity is deleted through the C API
   subroutine finalize_cavity_drop(self)
      type(cavity_type_drop), intent(inout) :: self

      !> Clean up timer
      ! TODO: optionally write timer summary (self%timer%write) before delete
      call self%timer%delete()

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
      if (allocated(self%v_i1_rA)) deallocate (self%v_i1_rA)
      if (allocated(self%xi1_rA)) deallocate (self%xi1_rA)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      if (allocated(self%cpjac_scal0)) deallocate (self%cpjac_scal0)
      if (allocated(self%A_tot1_rA)) deallocate (self%A_tot1_rA)
      if (allocated(self%asph1_rA)) deallocate (self%asph1_rA)
      if (allocated(self%V_tot1_rA)) deallocate (self%V_tot1_rA)
      if (allocated(self%vsph1_rA)) deallocate (self%vsph1_rA)

      ! Deallocate CPCM data
      if (allocated(self%cpcm_q)) deallocate (self%cpcm_q)
      if (allocated(self%cpcm_pot)) deallocate (self%cpcm_pot)
      if (allocated(self%cpcm_source_charges)) deallocate (self%cpcm_source_charges)
      if (allocated(self%cpcm_gradient)) deallocate (self%cpcm_gradient)

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
      if (allocated(self%qat)) deallocate (self%qat)
      if (allocated(self%aat)) deallocate (self%aat)
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%a)) deallocate (self%a)
      if (allocated(self%owner)) deallocate (self%owner)
      if (allocated(self%converged)) deallocate (self%converged)

   end subroutine finalize_cavity_drop

end module moist_cavity_drop
