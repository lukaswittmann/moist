!> Main DROP (Discretization via Reference-Onto surface Projection) implementation
module moist_cavity_drop
   use, intrinsic :: iso_fortran_env, only: output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_io_structure, only: structure_type
   use mctc_env, only: error_type, fatal_error, wp
   use moist_output_citations, only: citation
   use moist_math_linalg, only: setup_tangent_frame
   use moist_math_grid_angular_lebedev, only: lebedev_order_from_num
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, new_lebedev_grid
   use moist_cavity_type, only: cavity_type, list_cavity_fields_base, update_cavity_base, build_cavity_base, &
      & cavity_declare_pass
   use moist_channels_coupling, only: coupling_type, coupling_view_type, coupling_make_view, &
      & coupling_close_view
   use moist_channels_response, only: response_type, density_response_type, response_accumulate
   use moist_cavity_adjoint, only: cavity_adjoint_type
   use moist_channels_fields, only: field_query_type
   use moist_context, only: moist_context_type
   use moist_radius_type, only: radius_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   use moist_cavity_drop_gaussian, only: moist_cavity_drop_iswig, new_iswig
   use moist_math_adjacency, only: adjacency_list_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_sses, only: moist_cavity_drop_lsf_sses_type
   use moist_cavity_drop_lsf_isodensity, only: &
      moist_cavity_drop_lsf_isodensity_internal_type, &
      moist_cavity_drop_lsf_isodensity_callback_type
   use moist_cavity_drop_anchor, only: drop_anchor_type
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, drop_surface_weights_type
   use moist_math_smoothing, only: smoothing_kernel_wendland_type

   use moist_utils_timer, only: cat_setup, cat_solve, cat_properties, cat_gradient

   implicit none(type, external)
   private

   public :: cavity_type_drop
   public :: new_cavity_drop

   !> DROP cavity type
   type, extends(cavity_type) :: cavity_type_drop

      !> DROP parameters (blend_k, blend_1b, blend_2b, blend_3b, etc.)
      type(moist_cavity_drop_parameters_type) :: param

      !> Level set function model
      class(moist_cavity_drop_lsf_type), allocatable :: lsf_model

      !* ----------------------------- Atomic sphere data ----------------------------- *!

      !> Grid point neighbour list for density computation (CSR format)
      type(adjacency_list_type) :: grid_adj_list

      !> Unique numbering, `anchor%id + nsph*num_leb*(branch - 1)` (ngrid)
      integer, allocatable :: numbering(:)

      !* ------------------------- Anchor/row layout (CSR) ---------------------------- *!

      !> Unprojected anchor points, indexed by anchor (`anchor%ngrid`)
      type(drop_anchor_type) :: anchor

      !> Rows of anchor `i` are `offset(i)+1 : offset(i+1)` (anchor%ngrid + 1)
      !>
      !> - `anchor` arrays are sized `anchor%ngrid` and hold one value per
      !>   anchor; row arrays are sized `ngrid`, one per surface point
      !> - the base-type `owner` and `f` stay per row: they are the generic
      !>   cavity contract every model component reads
      integer, allocatable :: offset(:)
      !> Anchor of each row (ngrid)
      integer, allocatable :: anchor_of_row(:)

      !* ------------------------------ Integration grid ------------------------------ *!

      !> Final Lebedev quadrature weights (ngrid)
      real(wp), allocatable :: wleb(:)
      !> Raw Lebedev quadrature weights (3, nsph, ngrid)
      real(wp), allocatable :: wleb1_rA(:, :, :)

      !> Grid point density scaling
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

      !> Mean curvature
      real(wp), allocatable :: KM(:)
      !> Gaussian curvature
      real(wp), allocatable :: KG(:)

      !> Nuclear gradient of the first principal curvature (3, nsph, ngrid)
      !> Mean/Gaussian curvature gradients are derived from k1_rA/k2_rA
      !> downstream: dKM = (k1_rA + k2_rA)/2, dKG = k2*k1_rA + k1*k2_rA
      real(wp), allocatable :: k1_rA(:, :, :)
      !> Nuclear gradient of the second principal curvature (3, nsph, ngrid)
      real(wp), allocatable :: k2_rA(:, :, :)

      !* ------------------------------------ Areas ----------------------------------- *!

      !> Per-atom areas from the level set's atomic shares at each point (nsph)
      real(wp), allocatable :: asph_partition(:)
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

      !> Lebedev grid order index
      integer, allocatable :: oleb
      !> Cached Lebedev angular grid (3, num_leb)
      real(wp), allocatable :: ang_grid(:, :)
      !> Cached Lebedev solid-angle weights (num_leb), summing to 4*pi
      real(wp), allocatable :: ang_weight(:)

      !* ----------------------------- Switching functions ---------------------------- *!

      !> Critical level set switching function
      type(moist_math_smoothstep_sigmoid_type) :: f_crit
      !> Focal/branching point switching function
      type(moist_math_smoothstep_sigmoid_type) :: f_foc
      !> Lebedev weight pruning function (suppresses near-zero weights)
      type(moist_math_smoothstep_sigmoid_type) :: f_wleb

      !> Weight switching values per point (ngrid)
      real(wp), allocatable :: w_f0(:)
      !> Individual switching function component derivatives (3, nsph, ngrid)
      real(wp), allocatable :: w_f1_rA(:, :, :)

      !> iSwiG
      type(moist_cavity_drop_iswig) :: iswig

      !* --------------------------------- Projection --------------------------------- *!

      !> Branch weight model for handling degenerate branches
      type(branch_weight_type) :: branch_weight

      !> Canonical branch index of each row within its anchor (ngrid)
      integer, allocatable :: branch(:)
      !> Branch weight per grid point (ngrid) (1 = default/unbranched)
      real(wp), allocatable :: wbranch(:)
      !> Projection-objective value phi at each surviving grid point (ngrid)
      !> Used by the branch-weight softmax gradient to recover per-branch
      !> phi values once the main-loop KKT solve has produced dr/dR_A
      real(wp), allocatable :: phi0(:)

      real(wp), allocatable :: lambda0(:)
      !> Per-gridpoint derivatives of lambda w.r.t. nuclear coordinates (3, nsph, ngrid)
      real(wp), allocatable :: lambda1_rA(:, :, :)

   contains
      !> Enable optional geometry required by model components
      procedure :: configure => configure_cavity_drop
      !> Build the cavity for a new geometry
      procedure :: update => update_cavity_drop
      !> Whether the level set requires build-phase answers before the projection
      procedure :: defers_build => drop_defers_build
      !> Read the build-phase answers and project the deferred surface
      procedure :: build => build_cavity_drop
      !> Declare the build-phase requests of the level set
      procedure :: declare_coupling => declare_coupling_drop
      !> Compute area gradient w.r.t. nuclear coordinates
      procedure :: get_gradient => get_gradient_drop

      !> Compute CPCM properties
      procedure :: compute_gaussians

      !> Contract surface-coordinate weights to LSF adjoint weights
      procedure :: contract_surface_lsf_weights
      !> Map surface-coordinate weights into generic response channels
      procedure :: get_surface_response => get_surface_response_drop
      !> The DROP surface responds to the host field: host surface weights are
      !> consumed in the response phase, not only in the gradient
      procedure :: has_field_dependent_geometry => drop_field_dependent
      !> Contract surface-coordinate weights into the nuclear gradient
      procedure :: contract_surface_gradient => get_surface_gradient_drop

      !> Compute all needed cavity gradients
      procedure :: compute_gradient_drop

      !> Compute anchor-only nuclear derivatives (callback/isodensity LSF)
      procedure :: compute_anchor_gradient

      !> Build grid point neighbour list for density computation
      procedure :: setup_grid_adj_list

      !> Compute surface area from projected grid
      procedure :: compute_area_volume
      !> Partition the surface area among the atoms by the level set
      procedure :: compute_atom_area

      !> Cache Lebedev grid for current num_leb
      procedure :: ensure_lebedev_cache
      !> Project all grid points onto LSF surface
      procedure :: project_all_points
      !> Compute mean and Gaussian curvatures
      procedure :: compute_curvature
      !> Compute local grid point densities
      procedure :: compute_grid_point_density
      !> Run diagnostic checks on grid point health
      procedure :: analyze_cavity
      !> Write grid to CSV file for analysis
      procedure :: write_csv_debug => write_cavity_csv_debug
      !> Declare the readable results a DROP cavity holds
      procedure :: list_fields => list_cavity_fields_drop
      !> Construction settings
      procedure :: parameters => drop_parameters
      !> Print the cavity kind, the level set and the grouped settings
      procedure :: print_parameters => print_parameters_drop

      !> Finalizer
      final :: finalize_cavity_drop

   end type cavity_type_drop

   interface

      !* ============================================================================== *!
      !*              Internal DROP routines (should not be used outside)               *!
      !* ============================================================================== *!

      !> [projection.f90] Project all current grid points onto the LSF surface
      module subroutine project_all_points(self, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine project_all_points

      !> [properties.f90] Compute local grid point density values
      module subroutine compute_grid_point_density(self, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_grid_point_density

      !> [properties.f90] Check cavity diagnostics
      module subroutine analyze_cavity(self, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine analyze_cavity

      !> [properties.f90] Build grid point neighbour list for density computation
      module subroutine setup_grid_adj_list(self, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine setup_grid_adj_list

      !> [properties.f90] Compute total/atomic area and volume contributions
      module subroutine compute_area_volume(self, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_area_volume

      !> [properties.f90] Partition the surface area among the atoms by the level set
      module subroutine compute_atom_area(self, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_atom_area

      !> [properties.f90] Compute mean and Gaussian curvature on the surface grid
      module subroutine compute_curvature(self, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_curvature

      !> [properties.f90] Compute CPCM Gaussian width parameters for grid points
      module subroutine compute_gaussians(self, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine compute_gaussians

      !> [deriv/forward.f90] Compute first nuclear derivatives for DROP quantities
      module subroutine compute_gradient_drop(self, error, anchor_only)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
         !> Restrict each grid point's active atom to its owner (anchor motion
         !> only); used for callback/isodensity LSFs whose field nuclear
         !> derivatives are identically zero
         logical, intent(in), optional :: anchor_only
      end subroutine compute_gradient_drop

      !> [deriv/forward.f90] Compute anchor-only nuclear derivatives (owner motion,
      !> frozen field), a thin wrapper over
      !> `compute_gradient_drop(anchor_only=.true.)`
      module subroutine compute_anchor_gradient(self, error)
         implicit none(type, external)
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
         implicit none(type, external)
         class(cavity_type_drop), intent(in) :: self
         type(cavity_adjoint_type), intent(in) :: acc
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
         implicit none(type, external)
         class(cavity_type_drop), intent(in) :: self
         type(cavity_adjoint_type), intent(in) :: acc
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
         implicit none(type, external)
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
         implicit none(type, external)
         class(cavity_type_drop), intent(in) :: self
         type(cavity_adjoint_type), intent(in) :: acc
         real(wp), intent(out) :: w_lsf0(:)
         real(wp), intent(out) :: w_lsf1(:, :)
         real(wp), intent(out) :: w_lsf2(:, :, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine contract_surface_lsf_weights

      !> [deriv/potential.f90] Map accumulated surface adjoints into the generic response
      module subroutine get_surface_response_drop(self, acc, response, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(inout) :: self
         type(cavity_adjoint_type), intent(in) :: acc
         type(response_type), intent(inout) :: response
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_surface_response_drop

      !> [deriv/nuclear.f90] Contract accumulated surface adjoints into dE/dR_A
      !>
      !> @param[in]    self     DROP cavity instance
      !> @param[in]    acc      Accumulated surface-observable adjoints
      !> @param[inout] gradient Nuclear-gradient accumulator (3, nsph)
      !> @param[out]   error    Error object
      module subroutine get_surface_gradient_drop(self, acc, gradient, error)
         implicit none(type, external)
         class(cavity_type_drop), intent(in) :: self
         type(cavity_adjoint_type), intent(in) :: acc
         real(wp), intent(inout) :: gradient(:, :)
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_surface_gradient_drop

   end interface

contains

   !> Enable optional geometry requested by a model component
   !>
   !> @param[in,out] self DROP cavity
   !> @param[in] normal Request normals
   !> @param[in] curvature Request principal curvatures
   !> @param[out] error Configuration error
   subroutine configure_cavity_drop(self, normal, curvature, error)
      !> DROP cavity
      class(cavity_type_drop), intent(inout) :: self
      !> Request normals
      logical, intent(in), optional :: normal
      !> Request principal curvatures
      logical, intent(in), optional :: curvature
      !> Configuration error
      type(error_type), allocatable, intent(out) :: error

      if (present(normal)) self%param%do_normal = self%param%do_normal .or. normal
      if (present(curvature)) self%param%do_curvature = self%param%do_curvature .or. curvature
   end subroutine configure_cavity_drop

   !> The DROP surface geometry depends on the host field
   !>
   !> `get_surface_response_drop` maps the accumulated surface weights to a
   !> level-set response, so the host surface weights are consumed in the
   !> response phase and must be declared for it
   !>
   !> @param[in] self Cavity instance
   function drop_field_dependent(self) result(field_dependent)
      !> Cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Whether the live level set has a density adjoint
      logical :: field_dependent
      real(wp) :: factor
      field_dependent = .false.
      if (allocated(self%lsf_model)) field_dependent = self%lsf_model%density_adjoint_factor(factor)

   end function drop_field_dependent

   !> Construct from parameter values; omission uses compiled defaults
   !>
   !> @param[inout] self Object to initialize
   !> @param[in] radius_model Atomic radius model to copy
   !> @param[in] lsf_model Level set function to copy
   !> @param[out] error Construction error
   !> @param[in] param Configuration copied by value
   !> @param[in] ctx Borrowed run context; omitted, a model supplies its own
   subroutine new_cavity_drop(self, radius_model, lsf_model, error, param, ctx)
      !> Cavity to initialize
      type(cavity_type_drop), intent(inout) :: self
      !> Atomic radius model to copy
      class(radius_type), intent(in) :: radius_model
      !> Level set function to copy
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_model
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Configuration; omitted means compiled defaults
      type(moist_cavity_drop_parameters_type), intent(in), optional :: param
      !> Borrowed run context; omitted, a model supplies its own
      type(moist_context_type), intent(in), target, optional :: ctx

      nullify (self%ctx)
      if (present(ctx)) self%ctx => ctx

      call self%param%init_defaults()
      if (present(param)) self%param = param
      call self%param%compute_derived(error)
      if (allocated(error)) return

      !> Radius model setup
      if (allocated(self%radius_model)) deallocate (self%radius_model)
      allocate (self%radius_model, source=radius_model)

      !> LSF model setup
      if (allocated(self%lsf_model)) deallocate (self%lsf_model)
      allocate (self%lsf_model, source=lsf_model)
      !> Label by the level-set model, e.g. `SvdW-DROP`
      if (allocated(self%lsf_model%name)) then
         self%label = self%lsf_model%name//"-DROP"
      else
         self%label = "DROP"
      end if
      !> Push the cavity-derived screening threshold and cell-list cutoff into the LSF
      self%lsf_model%screening_threshold = self%param%screening_threshold
      call self%lsf_model%set_full_scan_below(self%param%cell_grid_full_scan_below)

      !> Set up weight switching function
      call new_smoothstep_sigmoid(self%f_crit, self%param%w_0ls_from, self%param%w_0ls_to, &
                                 p_hi=self%param%w_0ls_p, a_hi=self%param%w_0ls_a, p_lo=self%param%w_0ls_p, &
                                 a_lo=self%param%w_0ls_a)

      call new_smoothstep_sigmoid(self%f_foc, self%param%w_0tra_from, self%param%w_0tra_to, &
                                 p_hi=self%param%w_0ls_p, a_hi=self%param%w_0ls_a, p_lo=self%param%w_0ls_p, &
                                 a_lo=self%param%w_0ls_a)

      !> Set up Lebedev weight switching function (optional)
      if (self%param%wleb_prune_level > 0) then
         call new_smoothstep_sigmoid(self%f_wleb, &
                                    self%param%wleb_prune_from, self%param%wleb_prune_to)
      end if

      !> Set up hard-sphere gaussians function (iSwiG)
      call new_iswig(self%iswig, self%param%iswig_xi_born)

      !> Set up branch weight model
      call self%branch_weight%init(self%param%branch_weight_s)

      ! A constructed cavity has no structure, no surface and no coupling yet
      if (allocated(self%mol)) deallocate (self%mol)
      self%updated = .false.
      call self%clear_surface()
      call self%clear_couplings()

   end subroutine new_cavity_drop

   !> Record DROP and its level set on the cavity's context; called by
   !> `update`, which also covers a cavity adopted by a model
   subroutine cite_cavity_drop(self)
      !> Cavity instance
      class(cavity_type_drop), intent(inout) :: self

      if (.not. associated(self%ctx)) return
      call self%ctx%cite(citation%drop)
      select type (m => self%lsf_model)
      type is (moist_cavity_drop_lsf_cfc_type)
         call self%ctx%cite(citation%cfc)
      type is (moist_cavity_drop_lsf_isodensity_internal_type)
         call self%ctx%cite(citation%isodensity)
      type is (moist_cavity_drop_lsf_isodensity_callback_type)
         call self%ctx%cite(citation%isodensity)
      end select

   end subroutine cite_cavity_drop

   !> Construction settings of the cavity
   !>
   !> @param[in] self DROP cavity instance
   function drop_parameters(self) result(param)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in), target :: self
      !> Settings, valid while the cavity is
      class(moist_model_parameters_type), pointer :: param

      param => self%param

   end function drop_parameters

   !> Print the cavity section with the grouped level-set and DROP printers
   !>
   !> - headed by the label, e.g. `SvdW-DROP`
   !> - the level-set settings come first, then the DROP settings
   !> - the optional-property flags are part of the DROP settings
   !>
   !> @param[in] self DROP cavity instance
   !> @param[in] unit Output unit; defaults to the run context's unit
   !> @param[in] level Verbosity level; below 3 only the normal settings print,
   !>                  nothing at 0 (default: the context's
   !>                  verbosity, every setting without a context)
   subroutine print_parameters_drop(self, unit, level)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in), target :: self
      !> Output unit
      integer, intent(in), optional :: unit
      !> Verbosity level of the printout
      integer, intent(in), optional :: level

      !> Header printer
      type(prettyprinter) :: pp
      !> Level-set kind shown in the header
      character(len=:), allocatable :: lsf_kind
      !> Effective output unit and verbosity level
      integer :: iu, lvl

      iu = output_unit
      if (associated(self%ctx)) iu = self%ctx%unit
      if (present(unit)) iu = unit
      lvl = huge(1)
      if (associated(self%ctx)) lvl = self%ctx%verbosity()
      if (present(level)) lvl = level
      ! Silent at 0; from 1 an explicit request shows at least the normal settings
      if (lvl < 1) return
      lvl = max(2, lvl)
      lsf_kind = "none"
      if (allocated(self%lsf_model)) then
         lsf_kind = "unnamed"
         if (allocated(self%lsf_model%name)) lsf_kind = self%lsf_model%name
         ! The two isodensity models share a name; tell the density sources apart
         select type (m => self%lsf_model)
         type is (moist_cavity_drop_lsf_isodensity_internal_type)
            lsf_kind = lsf_kind//" (internal)"
         type is (moist_cavity_drop_lsf_isodensity_callback_type)
            lsf_kind = lsf_kind//" (callback)"
         end select
      end if

      pp = new_prettyprinter(unit=iu)
      if (allocated(self%label)) then
         call pp%push("Cavity ("//self%label//"):")
      else
         call pp%push("Cavity:")
      end if
      call pp%kv("Level set", lsf_kind)
      call pp%pop()
      call pp%blank()

      ! The level set the cavity is built on, then the DROP settings
      if (allocated(self%lsf_model)) then
         select type (m => self%lsf_model)
         type is (moist_cavity_drop_lsf_svdw_type)
            call m%param%print(unit=iu)
         type is (moist_cavity_drop_lsf_cfc_type)
            call m%param%print(unit=iu)
         type is (moist_cavity_drop_lsf_sses_type)
            call m%param%print(unit=iu)
         type is (moist_cavity_drop_lsf_isodensity_internal_type)
            call m%param%print(unit=iu, level=lvl)
         type is (moist_cavity_drop_lsf_isodensity_callback_type)
            call m%param%print(unit=iu, level=lvl)
         end select
      end if
      call self%param%print(unit=iu, level=lvl)

   end subroutine print_parameters_drop

   !* ================================================================================= *!
   !*                          Update Cavity (construct cavity)                         *!
   !* ================================================================================= *!

   !> Take a new geometry: radii and sphere centers, then the surface unless
   !> the level set waits for build-phase answers
   !>
   !> @param[inout] self     Cavity instance
   !> @param[in]    mol      Molecular structure
   !> @param[out]   error    Error handling
   subroutine update_cavity_drop(self, mol, error)
      class(cavity_type_drop), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error

      call update_cavity_base(self, mol, error)
      if (allocated(error)) return
      call self%require_context(error)
      if (allocated(error)) return
      call cite_cavity_drop(self)

      !> Set number of spheres
      self%nsph = mol%nat
      if (allocated(self%k1)) deallocate (self%k1)
      if (allocated(self%k2)) deallocate (self%k2)
      if (allocated(self%KM)) deallocate (self%KM)
      if (allocated(self%KG)) deallocate (self%KG)
      if (allocated(self%asph_partition)) deallocate (self%asph_partition)

      call self%radius_model%update(mol, error)
      if (allocated(error)) return
      if (allocated(self%radii)) deallocate (self%radii)
      allocate (self%radii(self%nsph))
      self%radii = self%radius_model%f0

      !> Set centers of spheres
      if (allocated(self%sphxyz)) deallocate (self%sphxyz)
      allocate (self%sphxyz(3, self%nsph), source=mol%xyz)

      !> The projection waits for the density matrix or the nuclear charges
      if (.not. self%defers_build()) then
         call project_surface_drop(self, error)
         if (allocated(error)) return
      end if
      self%updated = .true.
      self%built = .not. self%defers_build()

   end subroutine update_cavity_drop

   !> Whether the level set requires build-phase answers before the projection
   !>
   !> @param[in] self Cavity instance
   function drop_defers_build(self) result(defers)
      !> Cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Whether `update` leaves the projection to `build`
      logical :: defers

      defers = .false.
      if (allocated(self%lsf_model)) defers = self%lsf_model%requires_build_inputs()

   end function drop_defers_build

   !> Read the build phase, then project a deferred surface
   !>
   !> - `build_cavity_base` first; a surface that `update` projected stays
   !> - otherwise built once the level set read its answers and the
   !>   projection succeeded, then the coupling records the new point count
   !>
   !> @param[inout] self     Updated cavity
   !> @param[inout] coupling Coupling staged by `prepare_build`
   !> @param[out]   error    Error handling
   subroutine build_cavity_drop(self, coupling, error)
      !> Cavity instance
      class(cavity_type_drop), intent(inout) :: self
      !> Coupling staged by `prepare_build`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Cavity-scoped read interface
      type(coupling_view_type) :: view

      call build_cavity_base(self, coupling, error)
      if (allocated(error)) return
      if (.not. self%defers_build()) return
      call coupling_make_view(coupling, 0, view)
      call self%lsf_model%read_build(view, self%mol, error)
      call coupling_close_view(coupling)
      if (allocated(error)) return
      call project_surface_drop(self, error)
      if (allocated(error)) return
      self%built = .true.
      ! The answers of the build phase survive the new point count
      call cavity_declare_pass(self, coupling, error)

   end subroutine build_cavity_drop

   !> Project the surface of the stored structure, radii and sphere centers
   !>
   !> @param[inout] self  Cavity instance
   !> @param[out]   error Error handling
   subroutine project_surface_drop(self, error)
      !> Cavity instance
      class(cavity_type_drop), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Timer stack depth at entry; error paths unwind back to it (below)
      integer :: d0
      !> Anchors dropped by the pre-filter, and those below the iSwiG cut
      integer :: nremoved, n_iswig_removed
      !> Pre-filter report
      type(prettyprinter) :: pp

      ! The previous surface and its derivatives go before the projection
      ! re-allocates them: the getters test `A_tot1_rA` for availability
      call self%clear_surface()
      ! DROP row arrays; with `clear_surface` every array the projection allocates
      if (allocated(self%offset)) deallocate (self%offset)
      if (allocated(self%anchor_of_row)) deallocate (self%anchor_of_row)
      if (allocated(self%numbering)) deallocate (self%numbering)
      if (allocated(self%lambda0)) deallocate (self%lambda0)
      if (allocated(self%branch)) deallocate (self%branch)
      if (allocated(self%wbranch)) deallocate (self%wbranch)
      if (allocated(self%phi0)) deallocate (self%phi0)
      if (allocated(self%rho)) deallocate (self%rho)
      if (allocated(self%r_iI0)) deallocate (self%r_iI0)
      if (allocated(self%cpjac_scal0)) deallocate (self%cpjac_scal0)
      if (allocated(self%w_f0)) deallocate (self%w_f0)
      if (allocated(self%wleb)) deallocate (self%wleb)
      if (allocated(self%wleb1_rA)) deallocate (self%wleb1_rA)
      if (allocated(self%cpjac_scal1_rA)) deallocate (self%cpjac_scal1_rA)
      if (allocated(self%rho1_rA)) deallocate (self%rho1_rA)
      if (allocated(self%r_iI1_rA)) deallocate (self%r_iI1_rA)
      if (allocated(self%normal1_rA)) deallocate (self%normal1_rA)
      if (allocated(self%k1_rA)) deallocate (self%k1_rA)
      if (allocated(self%k2_rA)) deallocate (self%k2_rA)
      if (allocated(self%a_i1_rA)) deallocate (self%a_i1_rA)
      if (allocated(self%asph1_rA)) deallocate (self%asph1_rA)
      if (allocated(self%A_tot1_rA)) deallocate (self%A_tot1_rA)
      if (allocated(self%vsph1_rA)) deallocate (self%vsph1_rA)
      if (allocated(self%V_tot1_rA)) deallocate (self%V_tot1_rA)
      if (allocated(self%w_f1_rA)) deallocate (self%w_f1_rA)
      if (allocated(self%lambda1_rA)) deallocate (self%lambda1_rA)
      if (allocated(self%asph_partition)) deallocate (self%asph_partition)
      ! An optional list of the previous surface must not outlive it
      call self%grid_adj_list%destroy()

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

      !> Anchor grid: the Lebedev points on every sphere
      call self%ctx%timer%start("Array setup")
      call self%anchor%fill(self%sphxyz, self%radii, self%ang_grid, self%ang_weight)
      call self%ctx%timer%stop("Array setup")

      !> iSwiG atom-atom neighbour list, sorted by distance for early exit
      call self%ctx%timer%start("Adj. lists")
      call self%iswig%update(self%mol, self%radii, wleb_max=maxval(self%anchor%wleb0), &
                             nthreads=self%ctx%get_num_threads())
      call self%ctx%timer%stop("Adj. lists")

      !> Compute switching function
      call self%ctx%timer%start("Switching func.")
      call self%anchor%compute_switching(self%iswig, self%ctx%get_num_threads())
      if (self%ctx%verbosity(4)) then
         n_iswig_removed = count(self%anchor%f0 < self%param%wleb_cut)
         pp = new_prettyprinter(unit=self%ctx%unit)
         call pp%blank()
         call pp%push("Pre-filter:")
         call pp%kv("Points before removal", self%anchor%ngrid)
         call pp%kv("iSwig removed", n_iswig_removed)
         call pp%kv("Points to project", self%anchor%ngrid - n_iswig_removed)
         call pp%pop()
      end if
      call self%ctx%timer%stop("Switching func.")

      !> Pre-filter points below switching cutoff (reduces projection workload)
      call self%ctx%timer%start("Pre-filter")
      call self%anchor%filter(self%param%wleb_cut, nremoved)
      if (self%ctx%verbosity(4)) then
         pp = new_prettyprinter(unit=self%ctx%unit)
         call pp%blank()
         call pp%kv("Prefilter removed", nremoved, "points")
      end if
      call self%ctx%timer%stop("Pre-filter")

      call self%ctx%timer%stop("Setup")

      !* ------------------------------ Projection phase ------------------------------ *!
      call self%ctx%timer%start("Projector", category=cat_solve)
      call self%project_all_points(error)
      if (.not. allocated(error)) then
         call self%lsf_model%check_surface(self%xyz(:, 1:self%ngrid), &
                                           self%f(1:self%ngrid) >= self%param%wleb_cut, error)
      end if
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("Projector")

      !* ---------------------------- Post processing phase --------------------------- *!
      call self%ctx%timer%start("Post processing", category=cat_setup)

      ! CP Jacobian scaling, switching, branch weights and the weight cut ran
      ! per anchor inside the projection; only rows that carry weight exist

      !> Setup grid point adjacency list for density computation [optional]
      if (self%param%do_grid_density) then
         call self%ctx%timer%start("Grid adj. list")
         call self%setup_grid_adj_list(error)
         if (allocated(error)) then
            call self%ctx%timer%unwind(d0)
            return
         end if
         call self%ctx%timer%stop("Grid adj. list")
      end if

      call self%ctx%timer%stop("Post processing")

      !* ------------------------------ Properties phase ------------------------------ *!
      call self%ctx%timer%start("Properties", category=cat_properties)

      !> Compute grid point densities [optional diagnostic]
      if (self%param%do_grid_density) then
         call self%ctx%timer%start("Grid density")
         call self%compute_grid_point_density(error)
         if (allocated(error)) then
            call self%ctx%timer%unwind(d0)
            return
         end if
         call self%ctx%timer%stop("Grid density")
      end if

      !> Compute curvatures [optional diagnostic]
      if (self%param%do_curvature) then
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

      !> Partition the area among the atoms [optional diagnostic]
      if (self%param%do_atomarea) then
         call self%ctx%timer%start("Atom areas")
         call self%compute_atom_area(error)
         if (allocated(error)) then
            call self%ctx%timer%unwind(d0)
            return
         end if
         call self%ctx%timer%stop("Atom areas")
      end if

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
      if (self%ctx%verbosity(3)) then
         call self%ctx%timer%start("Analysis", category=cat_properties)
         call self%analyze_cavity(error)
         if (allocated(error)) then
            call self%ctx%timer%stop("Analysis")
            return
         end if
         call self%ctx%timer%stop("Analysis")
      end if

   end subroutine project_surface_drop

   !> Declare the build-phase requests of the level set at the cavity scope
   !>
   !> @param[in]    self     Cavity instance
   !> @param[inout] coupling Coupling being declared
   !> @param[out]   error    Error handling
   subroutine declare_coupling_drop(self, coupling, error)
      !> Cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Coupling being declared
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (allocated(self%lsf_model)) call self%lsf_model%declare_coupling(coupling, error)

   end subroutine declare_coupling_drop

   !* ================================================================================= *!
   !*                                 First derivatives                                 *!
   !* ================================================================================= *!

   !> Compute and store all requested DROP nuclear derivatives
   !>
   !> @param[inout] self  Cavity instance receiving the derivatives
   !> @param[out]   error Error handling
   subroutine get_gradient_drop(self, error)
      class(cavity_type_drop), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Timer stack depth at entry; error paths unwind back to it (below)
      integer :: d0

      call self%require_context(error)
      if (allocated(error)) return
      call self%require_built(error)
      if (allocated(error)) return

      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("Gradients", category=cat_gradient)

      if (self%ctx%verbosity(3)) write (self%ctx%unit, "(a)") "[Info] Computing gradients ..."
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
   !>
   !> - Unsupported sizes and rules with negative weights are errors of the
   !>   Lebedev module
   !> - Caches the nodes and solid-angle weights (summing to 4*pi) of the
   !>   unit-sphere grid
   !>
   !> @param[in,out] self   Cavity whose angular cache is refreshed
   !> @param[out]    error  Set for an unsupported or rejected Lebedev size
   subroutine ensure_lebedev_cache(self, error)
      !> Cavity instance
      class(cavity_type_drop), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Lebedev order index of the requested size
      integer :: oleb
      !> Unit-sphere grid whose nodes and weights are moved into the cache
      type(moist_math_grid_angular_type) :: leb

      ! Map requested num_leb to Lebedev order index
      call lebedev_order_from_num(self%param%num_leb, oleb, error, positive_weights_only=.true.)
      if (allocated(error)) return

      if (allocated(self%ang_grid) &
          .and. allocated(self%ang_weight) &
          .and. allocated(self%oleb) &
          ) then
         if (self%oleb == oleb .and. size(self%ang_weight) == self%param%num_leb) return
      end if

      if (allocated(self%ang_grid)) deallocate (self%ang_grid)
      if (allocated(self%ang_weight)) deallocate (self%ang_weight)
      if (allocated(self%oleb)) deallocate (self%oleb)

      call new_lebedev_grid(leb, error, npts=self%param%num_leb, positive_weights_only=.true.)
      if (allocated(error)) return

      call move_alloc(leb%points, self%ang_grid)
      call move_alloc(leb%weights, self%ang_weight)

      allocate (self%oleb)
      self%oleb = oleb

   end subroutine ensure_lebedev_cache

   !* ================================================================================= *!
   !*                                       Debug                                       *!
   !* ================================================================================= *!

   !> Write grid to csv file (for debugging)
   !>
   !> Per-point properties that were not computed are written as NaN
   !>
   !> @param[in]  self      DROP cavity instance
   !> @param[in]  filename  Output CSV path
   !> @param[out] error     Error on missing grid or failed open
   subroutine write_cavity_csv_debug(self, filename, error)
      class(cavity_type_drop), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, i
      real(wp) :: nan
      real(wp) :: val_rho_scal, val_rho_grid, val_rho_grid_anchor, val_KM, val_KG, val_cpjac

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, "write_csv_debug: cavity grid not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_csv_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open CSV file for writing: "//trim(filename))
         return
      end if

      write (unit, "(a)") "ngrid,numbering,x,y,z,owner,area,switch_f,w_leb,rho_scal0,rad,r_iI0,rho,"// &
         "anch_x,anch_y,anch_z,n_x,n_y,n_z,rho_grid,rho_grid_anchor,KM,KG,cpjac_scal"

      nan = ieee_value(nan, ieee_quiet_nan)
      do i = 1, self%ngrid

         val_rho_scal = nan
         val_rho_grid = nan
         val_rho_grid_anchor = nan
         val_KM = nan
         val_KG = nan
         val_cpjac = nan

         if (allocated(self%rho_scal0)) val_rho_scal = self%rho_scal0(i)
         if (allocated(self%rho_grid)) val_rho_grid = self%rho_grid(i)
         if (allocated(self%rho_grid_anchor)) val_rho_grid_anchor = self%rho_grid_anchor(i)
         if (allocated(self%KM)) val_KM = self%KM(i)
         if (allocated(self%KG)) val_KG = self%KG(i)
         if (allocated(self%cpjac_scal0)) val_cpjac = self%cpjac_scal0(i)

         write (unit, '(i0,",",i0,22(",",g0))') &
            i, &
            self%numbering(i), &
            self%xyz(1, i), self%xyz(2, i), self%xyz(3, i), &
            self%owner(i), &
            self%a(i), &
            self%f(i), &
            self%wleb(i), val_rho_scal, &
            self%radii(self%owner(i)), &
            self%r_iI0(i), &
            self%rho(i), &
            self%anchor%xyz(:, self%anchor_of_row(i)), &
            self%normal0(1, i), self%normal0(2, i), self%normal0(3, i), &
            val_rho_grid, val_rho_grid_anchor, &
            val_KM, val_KG, val_cpjac
      end do
      close (unit)

      write (self%ctx%unit, "(a,1x,a)") "[Info] Wrote cavity grid to", trim(filename)

   end subroutine write_cavity_csv_debug

   !* ================================================================================= *!
   !*                                  Readable results                                 *!
   !* ================================================================================= *!

   !> Declare the results a DROP cavity holds, on top of the generic ones
   !>
   !> Arrays gated by the `do_*` property flags of `param` are simply not allocated when they
   !> were not asked for, so they are not declared and a caller asking for one
   !> gets an error rather than zeros
   !>
   !> @param[in]    self   DROP cavity instance
   !> @param[inout] query  Walker collecting or fetching the declarations
   subroutine list_cavity_fields_drop(self, query)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Walker collecting or fetching the declarations
      type(field_query_type), intent(inout) :: query
      !> Per-row views of anchor arrays, built for this request only
      integer, allocatable :: row_int(:)
      real(wp), allocatable :: row_real(:), row_real2(:, :)

      call list_cavity_fields_base(self, query)
      ! Everything below belongs to the surface of the latest build
      if (.not. self%built) return

      !> Grid layout and provenance
      call query%add_int_value("num_leb", "Lebedev points per sphere the grid was built with", &
         & self%param%num_leb)
      call query%add_int_value("nmax", "Anchor points surviving the pre-filter", self%anchor%ngrid)
      call query%add_int_scalar("oleb", "Index of the cached Lebedev grid order", self%oleb)

      !> Point identity. `numbering` is the stable id across rebuilds; the
      !> per-row views of anchor quantities are materialized for this request only
      call query%add_int("numbering", "Stable point id, anchor_id + nsph*num_leb*(branch - 1) (ngrid)", &
         & self%numbering)
      call query%add_int("offset", "Row offset per anchor, 0-based start; rows of anchor i are "// &
         & "offset(i)+1:offset(i+1) (nmax + 1)", self%offset)
      call query%add_int("anchor_of_row", "Anchor slot of each row, 0-based (ngrid)", self%anchor_of_row, &
         & zero_based=.true.)
      if (allocated(self%anchor%id)) then
         row_int = self%anchor%id(self%anchor_of_row(1:self%ngrid))
         call query%add_int("anchor_id", "Anchor group id shared by the branches of one point (ngrid)", row_int)
      end if
      call query%add_int("branch", "Branch index within the anchor group, 1 when unbranched (ngrid)", &
         & self%branch)
      if (allocated(self%offset)) then
         row_int = self%offset(self%anchor_of_row(1:self%ngrid) + 1) &
            & - self%offset(self%anchor_of_row(1:self%ngrid))
         call query%add_int("branch_count", "Surviving branches in the point's anchor group (ngrid)", row_int)
      end if
      call query%add_real("wbranch", "Branch weight, 1 when unbranched (ngrid)", self%wbranch)

      !> Quadrature and projection geometry
      call query%add_real("wleb", "Lebedev quadrature weight per point (ngrid)", self%wleb)
      if (allocated(self%anchor%wleb0)) then
         row_real = self%anchor%wleb0(self%anchor_of_row(1:self%ngrid))
         call query%add_real("anchor_wleb0", "Lebedev quadrature weight of the anchor point (ngrid)", row_real)
      end if
      if (allocated(self%anchor%xyz)) then
         row_real2 = self%anchor%xyz(:, self%anchor_of_row(1:self%ngrid))
         call query%add_real2("anchorxyz", "Unprojected anchor coordinates, bohr (3, ngrid)", row_real2)
      end if
      call query%add_real("rho", "Anchor-to-projected-point distance, bohr (ngrid)", self%rho)
      call query%add_real("r_iI0", "Owner sphere center to grid point distance, bohr (ngrid)", self%r_iI0)
      call query%add_real("phi0", "Projection objective value at the point (ngrid)", self%phi0)
      call query%add_real("lambda0", "Projection Lagrange multiplier at the point (ngrid)", self%lambda0)
      call query%add_real("cpjac_scal0", "Closest-point Jacobian scaling (ngrid)", self%cpjac_scal0)
      call query%add_real("rho_scal0", "Grid point density scaling (ngrid)", self%rho_scal0)

      !> Switching functions
      call query%add_real("w_f0", "Weight switching value per point (ngrid)", self%w_f0)
      if (allocated(self%anchor%f0)) then
         row_real = self%anchor%f0(self%anchor_of_row(1:self%ngrid))
         call query%add_real("iswig_f0", "iSwiG switching value per point (ngrid)", row_real)
      end if
      if (allocated(self%anchor%xi0)) then
         row_real = self%anchor%xi0(self%anchor_of_row(1:self%ngrid))
         call query%add_real("anchor_xi0", "Gaussian width at the anchor point (ngrid)", row_real)
      end if

      !> Diagnostics, present only when the matching property was requested
      call query%add_real("KM", "Mean curvature (ngrid)", self%KM)
      call query%add_real("KG", "Gaussian curvature (ngrid)", self%KG)
      call query%add_real("rho_grid", "Local grid point density (ngrid)", self%rho_grid)
      call query%add_real("rho_grid_anchor", "Local grid point density, hard sphere limit (ngrid)", &
         & self%rho_grid_anchor)

      !> Per-sphere results
      call query%add_real("vsph", "Volume per sphere, bohr**3 (nsph)", self%vsph)

   end subroutine list_cavity_fields_drop

   !* ================================================================================= *!
   !*                                     Finalizer                                     *!
   !* ================================================================================= *!

   !> Finalizer deallocating every allocatable component of `cavity_type_drop`
   !>
   !> - ensures proper cleanup when the cavity is deleted through the C API
   subroutine finalize_cavity_drop(self)
      type(cavity_type_drop), intent(inout) :: self

      !> The profiling timer lives on the borrowed run context (self%ctx) and is
      !> owned by the top-level caller, so it is not torn down here; callers
      !> wanting the accumulated tree call self%ctx%timer%write(...) beforehand

      ! Deallocate grid point data arrays
      if (allocated(self%xi0)) deallocate (self%xi0)
      if (allocated(self%f)) deallocate (self%f)
      if (allocated(self%wleb)) deallocate (self%wleb)
      if (allocated(self%rho)) deallocate (self%rho)
      if (allocated(self%r_iI0)) deallocate (self%r_iI0)
      if (allocated(self%normal0)) deallocate (self%normal0)
      if (allocated(self%numbering)) deallocate (self%numbering)
      if (allocated(self%v)) deallocate (self%v)
      if (allocated(self%vsph)) deallocate (self%vsph)

      ! Deallocate anchor arrays
      call self%anchor%clear()

      ! Deallocate grid settings/caching
      if (allocated(self%oleb)) deallocate (self%oleb)
      if (allocated(self%ang_grid)) deallocate (self%ang_grid)
      if (allocated(self%ang_weight)) deallocate (self%ang_weight)

      ! Deallocate gradient arrays
      if (allocated(self%lambda0)) deallocate (self%lambda0)
      if (allocated(self%branch)) deallocate (self%branch)
      if (allocated(self%offset)) deallocate (self%offset)
      if (allocated(self%anchor_of_row)) deallocate (self%anchor_of_row)
      if (allocated(self%wbranch)) deallocate (self%wbranch)
      if (allocated(self%phi0)) deallocate (self%phi0)
      if (allocated(self%xyz1_rA)) deallocate (self%xyz1_rA)
      if (allocated(self%rho1_rA)) deallocate (self%rho1_rA)
      if (allocated(self%lambda1_rA)) deallocate (self%lambda1_rA)
      if (allocated(self%r_iI1_rA)) deallocate (self%r_iI1_rA)
      if (allocated(self%normal1_rA)) deallocate (self%normal1_rA)
      if (allocated(self%a_i1_rA)) deallocate (self%a_i1_rA)
      if (allocated(self%v1_rA)) deallocate (self%v1_rA)
      if (allocated(self%xi1_rA)) deallocate (self%xi1_rA)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      if (allocated(self%wleb1_rA)) deallocate (self%wleb1_rA)
      if (allocated(self%w_f1_rA)) deallocate (self%w_f1_rA)
      if (allocated(self%cpjac_scal0)) deallocate (self%cpjac_scal0)
      if (allocated(self%cpjac_scal1_rA)) deallocate (self%cpjac_scal1_rA)
      if (allocated(self%k1_rA)) deallocate (self%k1_rA)
      if (allocated(self%k2_rA)) deallocate (self%k2_rA)
      if (allocated(self%A_tot1_rA)) deallocate (self%A_tot1_rA)
      if (allocated(self%asph1_rA)) deallocate (self%asph1_rA)
      if (allocated(self%V_tot1_rA)) deallocate (self%V_tot1_rA)
      if (allocated(self%vsph1_rA)) deallocate (self%vsph1_rA)

      ! Deallocate neighbour list
      call self%grid_adj_list%destroy()

      ! Deallocate inherited allocatable components from base cavity_type
      if (allocated(self%mol)) deallocate (self%mol)
      if (allocated(self%asph)) deallocate (self%asph)
      if (allocated(self%asph_partition)) deallocate (self%asph_partition)
      if (allocated(self%total_area)) deallocate (self%total_area)
      if (allocated(self%total_volume)) deallocate (self%total_volume)
      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%sphxyz)) deallocate (self%sphxyz)
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%a)) deallocate (self%a)
      if (allocated(self%owner)) deallocate (self%owner)

   end subroutine finalize_cavity_drop

end module moist_cavity_drop
