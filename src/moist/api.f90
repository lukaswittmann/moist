!> C bindings for model configuration, host coupling, and cavity diagnostics
module moist_api
   use moist_cavity_drop_lsf_isodensity_param, only: moist_cavity_drop_lsf_isodensity_param_type
   use moist_cavity_iswig, only: moist_cavity_iswig_parameters_type
   use moist_model_component_pcm_type, only: moist_pcm_parameters_type
   use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, &
      & c_f_pointer, c_funptr, c_int, c_int8_t, c_int64_t, c_loc, c_null_char, &
      & c_null_ptr, c_ptr, c_size_t, c_sizeof
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_structure, only: structure_type, new
   use moist_cavity_type, only: cavity_type
   use moist_model_type, only: solvation_model_type, solvation_model_component_type
   use moist_channels_coupling, only: coupling_type, coupling_request_type, &
      & gaussian_moment_request_type, current_request, answer_flat, output_name_len
   use moist_channels_response, only: response_type, response_item_type, &
      & response_name_len, potential_adjoint_response_type, density_response_type, &
      & gostshyp_amplitude_response_type, current_response_item
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, &
                                             assemble_pcm_amat_with_gradient, pcm_amat_surface_weights, &
                                             pcm_amat_nuclear_gradient
   use moist_model_component_pcm_electrostatics, only: &
      pcm_electrostatic_nuclear_gradient
   use moist_model_component_gostshyp, only: solvation_model_component_gostshyp, &
      & new_component_gostshyp
   use moist_model_component_pcm_cpcm, only: solvation_model_component_cpcm, &
      & new_component_cpcm
   use moist_model_component_pcm_cosmo, only: solvation_model_component_cosmo, &
      & new_component_cosmo
   use moist_model_component_pv, only: solvation_model_component_pv, new_component_pv
   use moist_model_general, only: solvation_model_general, new_model_general
   use moist_context, only: moist_context_type, new_context
   use moist_radii, only: radius_type, new_radii, radius_type_custom
   use moist_radii_custom, only: new_custom_radii_atoms, new_custom_radii_elements
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_lsf_cfc_param, only: moist_cavity_drop_lsf_cfc_param_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_isodensity_callback, only: &
      moist_cavity_drop_lsf_isodensity_callback_type
   use moist_cavity_drop_lsf_isodensity_internal, only: &
      moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_cavity_fields, only: cavity_field_query_type, cavity_field_max_rank, &
      & cavity_field_real, cavity_field_int, cavity_field_bool
   use moist_version, only: get_moist_version
   use moist_output_ascii, only: moist_banner_text
   implicit none(type, external)
   private

   character(len=*), parameter :: namespace = "moist_"
   integer, parameter :: api_max_cstr = 4096
   !> Longest cavity field name the API will scan for
   integer, parameter :: max_field_name_len = 64

   !> Frozen 1.0 drop layout; a future addition needs a new versioned type
   type, bind(C) :: api_drop_options_v1_0
      !> Caller allocation size in bytes
      integer(c_size_t) :: struct_size = 0_c_size_t
      !> Lebedev points per atom
      integer(c_int) :: nleb = 194_c_int
      !> Enable diagnostic checks
      logical(c_bool) :: debug = .false._c_bool
      !> Output detail level
      integer(c_int) :: verbosity = 0_c_int
      !> Enable fine cavity refinement
      logical(c_bool) :: do_fine = .false._c_bool
      !> Projection convergence tolerance
      real(c_double) :: tolerance = 1.0e-10_c_double
      !> Maximum projection iterations
      integer(c_int) :: proj_maxiter = 150_c_int
      !> Projection refinement level
      integer(c_int) :: proj_level = 3_c_int
      !> Branch-selection smoothing weight
      real(c_double) :: branch_weight_s = 0.0025_c_double
      !> Grid density kernel length in Bohr
      real(c_double) :: rho_grid_h = 1.0_c_double
      !> Lebedev pruning level
      integer(c_int) :: wleb_prune_level = 0_c_int
      !> Reserved ABI storage; initialized to zero and ignored on input
      integer(c_int) :: reserved0 = 0_c_int
   end type api_drop_options_v1_0

   !> Minimum valid caller size, independent of future library extensions
   integer(c_size_t), parameter :: drop_options_min_size = c_sizeof(api_drop_options_v1_0())

   !> Frozen 1.0 iswig layout; a future addition needs a new versioned type
   type, bind(C) :: api_iswig_options_v1_0
      !> Caller allocation size in bytes
      integer(c_size_t) :: struct_size = 0_c_size_t
      !> Lebedev points per atom
      integer(c_int) :: nleb = 110_c_int
      !> Enable diagnostic checks
      logical(c_bool) :: debug = .false._c_bool
      !> Output detail level
      integer(c_int) :: verbosity = 0_c_int
      !> Surface-weight cutoff
      real(c_double) :: cut_a = 0.0_c_double
      !> Switching-factor cutoff
      real(c_double) :: cut_f = 1.0e-10_c_double
   end type api_iswig_options_v1_0

   !> Minimum valid caller size, independent of future library extensions
   integer(c_size_t), parameter :: iswig_options_min_size = c_sizeof(api_iswig_options_v1_0())

   !> Frozen 1.0 svdw layout; a future addition needs a new versioned type
   type, bind(C) :: api_svdw_options_v1_0
      !> Caller allocation size in bytes
      integer(c_size_t) :: struct_size = 0_c_size_t
      !> Soft minimum sharpness
      real(c_double) :: blend_k = 5.5_c_double
      !> One-body blending weight
      real(c_double) :: blend_1b = 1.0_c_double
      !> Two-body blending weight
      real(c_double) :: blend_2b = 0.0_c_double
      !> Three-body blending weight
      real(c_double) :: blend_3b = 3.0_c_double
   end type api_svdw_options_v1_0

   !> Minimum valid caller size, independent of future library extensions
   integer(c_size_t), parameter :: svdw_options_min_size = c_sizeof(api_svdw_options_v1_0())

   !> Frozen 1.0 cfc layout; a future addition needs a new versioned type
   type, bind(C) :: api_cfc_options_v1_0
      !> Caller allocation size in bytes
      integer(c_size_t) :: struct_size = 0_c_size_t
      !> Atomic-term exponent
      real(c_double) :: a1 = -15.0_c_double
      !> Pair-term exponent
      real(c_double) :: a2 = -9.0_c_double
      !> Pair-term coupling coefficient
      real(c_double) :: c = 5.0_c_double
      !> Pair-term polynomial power
      integer(c_int) :: m = 4_c_int
      !> Reserved ABI storage; initialized to zero and ignored on input
      integer(c_int) :: reserved0 = 0_c_int
   end type api_cfc_options_v1_0

   !> Minimum valid caller size, independent of future library extensions
   integer(c_size_t), parameter :: cfc_options_min_size = c_sizeof(api_cfc_options_v1_0())

   !> Frozen 1.0 isodensity layout; a future addition needs a new versioned type
   type, bind(C) :: api_isodensity_options_v1_0
      !> Caller allocation size in bytes
      integer(c_size_t) :: struct_size = 0_c_size_t
      !> Electron-density isovalue in atomic units
      real(c_double) :: rho_iso = 1.0e-3_c_double
      !> Positive level-set scaling factor
      real(c_double) :: scale = 1000.0_c_double
   end type api_isodensity_options_v1_0

   !> Minimum valid caller size, independent of future library extensions
   integer(c_size_t), parameter :: isodensity_options_min_size = c_sizeof(api_isodensity_options_v1_0())

   !> Frozen 1.0 model layout; a future addition needs a new versioned type
   type, bind(C) :: api_model_options_v1_0
      !> Caller allocation size in bytes
      integer(c_size_t) :: struct_size = 0_c_size_t
      !> Enable diagnostic checks
      logical(c_bool) :: debug = .false._c_bool
      !> Output detail level
      integer(c_int) :: verbosity = 0_c_int
   end type api_model_options_v1_0

   !> Minimum valid caller size, independent of future library extensions
   integer(c_size_t), parameter :: model_options_min_size = c_sizeof(api_model_options_v1_0())

   !> Frozen 1.0 pcm layout; a future addition needs a new versioned type
   type, bind(C) :: api_pcm_options_v1_0
      !> Caller allocation size in bytes
      integer(c_size_t) :: struct_size = 0_c_size_t
      !> PCM linear solver selector
      integer(c_int) :: solver = 3_c_int
      !> Iteration cap; only read by the iterative solver
      integer(c_int) :: solver_maxiter = 50_c_int
      !> Residual threshold; only read by the iterative solver
      real(c_double) :: solver_tol = 1.0e-10_c_double
   end type api_pcm_options_v1_0

   !> Minimum valid caller size, independent of future library extensions
   integer(c_size_t), parameter :: pcm_options_min_size = c_sizeof(api_pcm_options_v1_0())

   !> Owning C handle for an independently configured level-set function
   type :: vp_lsf
      !> Concrete level set, copied by the cavity constructor
      class(moist_cavity_drop_lsf_type), allocatable :: ptr
   end type vp_lsf

   type :: vp_error
      type(error_type), allocatable :: ptr
   end type vp_error

   type :: vp_structure
      type(structure_type) :: ptr
   end type vp_structure

   type :: vp_cavity
      !> Run context owned by this handle; the cavity borrows a pointer to it
      type(moist_context_type) :: ctx
      class(cavity_type), pointer :: ptr => null()
      logical :: owned = .true.
   end type vp_cavity

   type :: vp_radii
      class(radius_type), allocatable :: ptr
   end type vp_radii

   type :: vp_model
      !> Run context owned by this handle; the model borrows a pointer to it
      !>
      !> - built by a model constructor with `new_context`, then passed to the
      !>   concrete model and on to its components, as the cavity handles below
      !> - torn down in `delete_solvation_model_api`
      type(moist_context_type) :: ctx
      class(solvation_model_type), allocatable :: ptr
   end type vp_model

   type :: vp_component
      !> Run context owned by this handle until the component is copied into a
      !> model; `solvation_model_general%add_component` re-points the copy at the
      !> model context so the copy stays valid after this handle is deleted
      type(moist_context_type) :: ctx
      !> Concrete component owned by this opaque handle
      class(solvation_model_component_type), allocatable :: ptr
   end type vp_component

   !> Borrowed host handle for a coupling owned by its general model
   !>
   !> - the model must outlive this wrapper
   !> - deletion releases its collection
   type :: vp_coupling
      !> Coupling owned by the parent model
      type(coupling_type), pointer :: ptr => null()
      !> Parent model whose registry owns the collection
      type(solvation_model_general), pointer :: owner => null()
   end type vp_coupling

   !> Response handle (`moist_response`), filled by the `get_*` reads
   !>
   !> Every valid read clears it and returns the complete host part of its
   !> phase; the host walks the items and copies the arrays of the current one
   !> out by name
   type :: vp_response
      !> Response list owned by this handle
      type(response_type) :: ptr
   end type vp_response

   !> Successful error status
   integer(c_int), parameter :: api_success = 0_c_int
   !> Error handle contains a diagnostic
   integer(c_int), parameter :: api_failure = 1_c_int
   !> Required error handle is absent
   integer(c_int), parameter :: api_invalid_error = 2_c_int

   public :: vp_error, vp_structure, vp_cavity, vp_radii, vp_model, vp_component
   public :: vp_coupling, vp_response
   public :: get_version_api, get_version_string_api
   public :: new_error_api, check_error_api, get_error_api, delete_error_api
   public :: new_structure_api, delete_structure_api, update_structure_api
   public :: new_cpcm_radii_api, new_smd_radii_api, new_d3_radii_api, new_cosmo_radii_api, new_bondi_radii_api
   public :: new_custom_radii_api, set_custom_radii_atoms_api, set_custom_radii_elements_api
   public :: delete_radii_api
   ! Solvation model API
   public :: update_solvation_model_api
   public :: get_solvation_model_cavity_api
   public :: delete_solvation_model_api
   ! General solvation model and its components
   public :: new_cpcm_component_api, new_cosmo_component_api
   public :: new_pv_component_api, new_gostshyp_component_api
   public :: delete_solvation_component_api
   public :: new_general_solvation_model_api, general_model_add_component_api
   ! Host coupling protocol (C mirror of the Fortran request/response exchange)
   public :: new_coupling_api, delete_coupling_api
   public :: new_response_api, delete_response_api
   public :: general_model_prepare_energy_api, general_model_prepare_response_api
   public :: general_model_prepare_gradient_api
   public :: general_model_get_energy_api, general_model_get_response_api
   public :: general_model_get_gradient_api
   public :: next_coupling_request_api, coupling_request_name_api
   public :: coupling_request_missing_api
   public :: coupling_answer_api, coupling_get_gaussian_moment_width_api
   public :: next_response_item_api, response_item_name_api, response_get_api
   ! Type-specific constructors

   ! Generic cavity operations
   public :: update_cavity_api
   public :: get_cavity_sizes_api
   public :: get_cavity_results_api
   public :: delete_cavity_api
   ! Named cavity result fields
   public :: get_cavity_field_count_api
   public :: get_cavity_field_info_api
   public :: get_cavity_field_about_api
   public :: get_cavity_field_real_api
   public :: get_cavity_field_int_api
   public :: get_cavity_field_bool_api
   ! Legacy DROP API (deprecated - use the generic cavity operations above)

   ! Type-specific getters
   public :: get_drop_cavity_tolerance_api
   public :: get_isodensity_cart_layout_api
   public :: set_isodensity_density_api
   public :: assemble_amat_api
   public :: get_cavity_gaussian_api
   ! Gradient API
   public :: compute_cavity_gradient_api
   public :: compute_anchor_gradient_api
   public :: get_anchor_gradient_api
   public :: get_cavity_gradient_api
   public :: get_amat_gradient_api
   public :: contract_amat1_q1q2_rA_api
   public :: contract_amat1_q1q2_surface_weights_api
   public :: contract_surface_lsf_weights_api, contract_surface_lsf_weights_extended_api
   public :: contract_pcm_nuclear_gradient_api

contains

   !> Obtain isodensity defaults from the native parameter type
   function default_isodensity_options() result(value)
      !> Interoperable options
      type(api_isodensity_options_v1_0) :: value
      !> Native compiled defaults
      type(moist_cavity_drop_lsf_isodensity_param_type) :: native
      value%rho_iso = native%rho_iso
      value%scale = native%scale
   end function default_isodensity_options

   !> Obtain pcm defaults from the native parameter type
   function default_pcm_options() result(value)
      !> Interoperable options
      type(api_pcm_options_v1_0) :: value
      !> Native compiled defaults
      type(moist_pcm_parameters_type) :: native
      value%solver = native%solver
      value%solver_maxiter = native%solver_maxiter
      value%solver_tol = native%solver_tol
   end function default_pcm_options

   !> Obtain cfc defaults from the native parameter type
   function default_cfc_options() result(value)
      !> Interoperable options
      type(api_cfc_options_v1_0) :: value
      !> Native compiled defaults
      type(moist_cavity_drop_lsf_cfc_param_type) :: native
      value%a1 = native%a1
      value%a2 = native%a2
      value%c = native%c
      value%m = native%m
   end function default_cfc_options

   !> Obtain svdw defaults from the native parameter type
   function default_svdw_options() result(value)
      !> Interoperable options
      type(api_svdw_options_v1_0) :: value
      !> Native compiled defaults
      type(moist_cavity_drop_lsf_svdw_param_type) :: native
      value%blend_k = native%blend_k
      value%blend_1b = native%blend_1b
      value%blend_2b = native%blend_2b
      value%blend_3b = native%blend_3b
   end function default_svdw_options

   !> Obtain iswig defaults from the native parameter type
   function default_iswig_options() result(value)
      !> Interoperable options
      type(api_iswig_options_v1_0) :: value
      !> Native compiled defaults
      type(moist_cavity_iswig_parameters_type) :: native
      value%nleb = native%num_leb
      value%cut_a = native%cut_a
      value%cut_f = native%cut_f
   end function default_iswig_options

   !> Obtain drop defaults from the native parameter type
   function default_drop_options() result(value)
      !> Interoperable options
      type(api_drop_options_v1_0) :: value
      !> Native compiled defaults
      type(moist_cavity_drop_parameters_type) :: native
      value%nleb = native%num_leb
      value%tolerance = native%tolerance
      value%proj_maxiter = native%proj_maxiter
      value%proj_level = native%proj_level
      value%branch_weight_s = native%branch_weight_s
      value%rho_grid_h = native%rho_grid_h
      value%wleb_prune_level = native%wleb_prune_level
      value%do_fine = logical(native%do_fine, c_bool)
   end function default_drop_options

   !* ================================================================================= *!
   !*                              Version and diagnostics                              *!
   !* ================================================================================= *!

   !> Copy banner text or query its length; printing belongs to the host
   subroutine get_banner_api(verror, style, buffer, capacity, length) bind(C, name="moist_get_banner")
      !> Required diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Banner style selector
      integer(c_int), value, intent(in) :: style
      !> Output buffer; NULL is allowed only with zero capacity
      character(kind=c_char), intent(inout), optional :: buffer(*)
      !> Buffer capacity in bytes, including the terminator
      integer(c_size_t), value, intent(in) :: capacity
      !> Required text length excluding the terminator; unchanged on error
      integer(c_size_t), intent(inout), optional :: length
      !> Decoded diagnostic handle
      type(vp_error), pointer :: error

      if (.not. valid_string_output(verror, buffer, capacity, length, "get_banner", error)) return
      if (style < 0 .or. style > 3) then
         call api_error(error%ptr, "get_banner", "Invalid banner style")
         return
      end if
      call copy_string_output(moist_banner_text(int(style)), buffer, capacity, length)
   end subroutine get_banner_api

   !> Validate a string query and terminate any writable output buffer
   function valid_string_output(verror, buffer, capacity, length, routine, error) result(valid)
      !> Required diagnostic handle
      type(c_ptr), intent(in) :: verror
      !> Caller output, cleared even when validation fails
      character(kind=c_char), intent(inout), optional :: buffer(*)
      !> Buffer capacity including the terminator
      integer(c_size_t), intent(in) :: capacity
      !> Required output pointer; its value is preserved on failure
      integer(c_size_t), intent(in), optional :: length
      !> Calling entry point for diagnostics
      character(len=*), intent(in) :: routine
      !> Decoded error handle when supplied
      type(vp_error), pointer, intent(out) :: error
      !> Whether copying is permitted
      logical :: valid

      valid = .false.
      nullify (error)
      if (present(buffer)) then
         if (capacity /= 0) buffer(1) = c_null_char
      end if
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(length)) then
         call api_error(error%ptr, routine, "Length output is required")
         return
      end if
      if (capacity < 0 .or. (capacity > 0 .and. .not. present(buffer))) then
         call api_error(error%ptr, routine, "Invalid buffer or capacity")
         return
      end if
      valid = .true.
   end function valid_string_output

   !> Copy a validated string query, truncating while reporting its full length
   subroutine copy_string_output(text, buffer, capacity, length)
      !> Full source text
      character(len=*), intent(in) :: text
      !> Optional destination for a size-only query
      character(kind=c_char), intent(inout), optional :: buffer(*)
      !> Caller allocation size in bytes
      integer(c_size_t), intent(in) :: capacity
      !> Full length, excluding the terminator
      integer(c_size_t), intent(out) :: length
      !> Number of characters copied and loop index
      integer(c_size_t) :: copied, i

      length = len(text, kind=c_size_t)
      if (capacity == 0 .or. .not. present(buffer)) return
      copied = min(length, capacity - 1_c_size_t)
      do i = 1, copied
         buffer(i) = text(i:i)
      end do
      buffer(copied + 1) = c_null_char
   end subroutine copy_string_output

   !* ================================================================================= *!
   !*                              Options and constructors                             *!
   !* ================================================================================= *!

   !> Initialize drop options within the caller's allocation
   subroutine init_drop_options_api(verror, options, bytes) bind(C, name="moist_init_drop_options")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Caller-owned options buffer
      type(c_ptr), value, intent(in) :: options
      !> Allocated byte count
      integer(c_size_t), value, intent(in) :: bytes
      !> Compiled defaults
      type(api_drop_options_v1_0), target :: defaults
      !> Default values; only named components are copied to the wire layout
      type(api_drop_options_v1_0) :: values
      !> Decoded error handle
      type(vp_error), pointer :: error
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      values = default_drop_options()
      call zero_options_storage(c_loc(defaults), c_sizeof(defaults))
      defaults%struct_size = bytes
      defaults%nleb = values%nleb
      defaults%debug = values%debug
      defaults%verbosity = values%verbosity
      defaults%do_fine = values%do_fine
      defaults%tolerance = values%tolerance
      defaults%proj_maxiter = values%proj_maxiter
      defaults%proj_level = values%proj_level
      defaults%branch_weight_s = values%branch_weight_s
      defaults%rho_grid_h = values%rho_grid_h
      defaults%wleb_prune_level = values%wleb_prune_level
      defaults%reserved0 = values%reserved0
      call copy_options(c_loc(defaults), options, c_sizeof(defaults), bytes, drop_options_min_size, error%ptr)
      call prefix_api_error(error%ptr, "init_drop_options")
   end subroutine init_drop_options_api

   !> Decode optional drop options into an independent value
   subroutine read_drop_options(options, value, error)
      !> Optional caller options; NULL selects defaults
      type(c_ptr), intent(in) :: options
      !> Default-initialized result
      type(api_drop_options_v1_0), intent(out), target :: value
      !> Diagnostic on invalid size
      type(error_type), allocatable, intent(out) :: error
      !> Caller allocation size, the common first field
      integer(c_size_t), pointer :: bytes
      value = default_drop_options()
      value%struct_size = c_sizeof(value)
      if (.not. c_associated(options)) return
      call c_f_pointer(options, bytes)
      call copy_options(options, c_loc(value), c_sizeof(value), bytes, drop_options_min_size, error)
   end subroutine read_drop_options

   !> Initialize iswig options within the caller's allocation
   subroutine init_iswig_options_api(verror, options, bytes) bind(C, name="moist_init_iswig_options")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Caller-owned options buffer
      type(c_ptr), value, intent(in) :: options
      !> Allocated byte count
      integer(c_size_t), value, intent(in) :: bytes
      !> Compiled defaults
      type(api_iswig_options_v1_0), target :: defaults
      !> Default values; only named components are copied to the wire layout
      type(api_iswig_options_v1_0) :: values
      !> Decoded error handle
      type(vp_error), pointer :: error
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      values = default_iswig_options()
      call zero_options_storage(c_loc(defaults), c_sizeof(defaults))
      defaults%struct_size = bytes
      defaults%nleb = values%nleb
      defaults%debug = values%debug
      defaults%verbosity = values%verbosity
      defaults%cut_a = values%cut_a
      defaults%cut_f = values%cut_f
      call copy_options(c_loc(defaults), options, c_sizeof(defaults), bytes, iswig_options_min_size, error%ptr)
      call prefix_api_error(error%ptr, "init_iswig_options")
   end subroutine init_iswig_options_api

   !> Decode optional iswig options into an independent value
   subroutine read_iswig_options(options, value, error)
      !> Optional caller options; NULL selects defaults
      type(c_ptr), intent(in) :: options
      !> Default-initialized result
      type(api_iswig_options_v1_0), intent(out), target :: value
      !> Diagnostic on invalid size
      type(error_type), allocatable, intent(out) :: error
      !> Caller allocation size, the common first field
      integer(c_size_t), pointer :: bytes
      value = default_iswig_options()
      value%struct_size = c_sizeof(value)
      if (.not. c_associated(options)) return
      call c_f_pointer(options, bytes)
      call copy_options(options, c_loc(value), c_sizeof(value), bytes, iswig_options_min_size, error)
   end subroutine read_iswig_options

   !> Initialize svdw options within the caller's allocation
   subroutine init_svdw_options_api(verror, options, bytes) bind(C, name="moist_init_svdw_options")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Caller-owned options buffer
      type(c_ptr), value, intent(in) :: options
      !> Allocated byte count
      integer(c_size_t), value, intent(in) :: bytes
      !> Compiled defaults
      type(api_svdw_options_v1_0), target :: defaults
      !> Default values; only named components are copied to the wire layout
      type(api_svdw_options_v1_0) :: values
      !> Decoded error handle
      type(vp_error), pointer :: error
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      values = default_svdw_options()
      call zero_options_storage(c_loc(defaults), c_sizeof(defaults))
      defaults%struct_size = bytes
      defaults%blend_k = values%blend_k
      defaults%blend_1b = values%blend_1b
      defaults%blend_2b = values%blend_2b
      defaults%blend_3b = values%blend_3b
      call copy_options(c_loc(defaults), options, c_sizeof(defaults), bytes, svdw_options_min_size, error%ptr)
      call prefix_api_error(error%ptr, "init_svdw_options")
   end subroutine init_svdw_options_api

   !> Decode optional svdw options into an independent value
   subroutine read_svdw_options(options, value, error)
      !> Optional caller options; NULL selects defaults
      type(c_ptr), intent(in) :: options
      !> Default-initialized result
      type(api_svdw_options_v1_0), intent(out), target :: value
      !> Diagnostic on invalid size
      type(error_type), allocatable, intent(out) :: error
      !> Caller allocation size, the common first field
      integer(c_size_t), pointer :: bytes
      value = default_svdw_options()
      value%struct_size = c_sizeof(value)
      if (.not. c_associated(options)) return
      call c_f_pointer(options, bytes)
      call copy_options(options, c_loc(value), c_sizeof(value), bytes, svdw_options_min_size, error)
   end subroutine read_svdw_options

   !> Initialize cfc options within the caller's allocation
   subroutine init_cfc_options_api(verror, options, bytes) bind(C, name="moist_init_cfc_options")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Caller-owned options buffer
      type(c_ptr), value, intent(in) :: options
      !> Allocated byte count
      integer(c_size_t), value, intent(in) :: bytes
      !> Compiled defaults
      type(api_cfc_options_v1_0), target :: defaults
      !> Default values; only named components are copied to the wire layout
      type(api_cfc_options_v1_0) :: values
      !> Decoded error handle
      type(vp_error), pointer :: error
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      values = default_cfc_options()
      call zero_options_storage(c_loc(defaults), c_sizeof(defaults))
      defaults%struct_size = bytes
      defaults%a1 = values%a1
      defaults%a2 = values%a2
      defaults%c = values%c
      defaults%m = values%m
      defaults%reserved0 = values%reserved0
      call copy_options(c_loc(defaults), options, c_sizeof(defaults), bytes, cfc_options_min_size, error%ptr)
      call prefix_api_error(error%ptr, "init_cfc_options")
   end subroutine init_cfc_options_api

   !> Decode optional cfc options into an independent value
   subroutine read_cfc_options(options, value, error)
      !> Optional caller options; NULL selects defaults
      type(c_ptr), intent(in) :: options
      !> Default-initialized result
      type(api_cfc_options_v1_0), intent(out), target :: value
      !> Diagnostic on invalid size
      type(error_type), allocatable, intent(out) :: error
      !> Caller allocation size, the common first field
      integer(c_size_t), pointer :: bytes
      value = default_cfc_options()
      value%struct_size = c_sizeof(value)
      if (.not. c_associated(options)) return
      call c_f_pointer(options, bytes)
      call copy_options(options, c_loc(value), c_sizeof(value), bytes, cfc_options_min_size, error)
   end subroutine read_cfc_options

   !> Initialize isodensity options within the caller's allocation
   subroutine init_isodensity_options_api(verror, options, bytes) bind(C, name="moist_init_isodensity_options")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Caller-owned options buffer
      type(c_ptr), value, intent(in) :: options
      !> Allocated byte count
      integer(c_size_t), value, intent(in) :: bytes
      !> Compiled defaults
      type(api_isodensity_options_v1_0), target :: defaults
      !> Default values; only named components are copied to the wire layout
      type(api_isodensity_options_v1_0) :: values
      !> Decoded error handle
      type(vp_error), pointer :: error
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      values = default_isodensity_options()
      call zero_options_storage(c_loc(defaults), c_sizeof(defaults))
      defaults%struct_size = bytes
      defaults%rho_iso = values%rho_iso
      defaults%scale = values%scale
      call copy_options(c_loc(defaults), options, c_sizeof(defaults), bytes, isodensity_options_min_size, error%ptr)
      call prefix_api_error(error%ptr, "init_isodensity_options")
   end subroutine init_isodensity_options_api

   !> Decode optional isodensity options into an independent value
   subroutine read_isodensity_options(options, value, error)
      !> Optional caller options; NULL selects defaults
      type(c_ptr), intent(in) :: options
      !> Default-initialized result
      type(api_isodensity_options_v1_0), intent(out), target :: value
      !> Diagnostic on invalid size
      type(error_type), allocatable, intent(out) :: error
      !> Caller allocation size, the common first field
      integer(c_size_t), pointer :: bytes
      value = default_isodensity_options()
      value%struct_size = c_sizeof(value)
      if (.not. c_associated(options)) return
      call c_f_pointer(options, bytes)
      call copy_options(options, c_loc(value), c_sizeof(value), bytes, isodensity_options_min_size, error)
   end subroutine read_isodensity_options

   !> Initialize model options within the caller's allocation
   subroutine init_model_options_api(verror, options, bytes) bind(C, name="moist_init_model_options")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Caller-owned options buffer
      type(c_ptr), value, intent(in) :: options
      !> Allocated byte count
      integer(c_size_t), value, intent(in) :: bytes
      !> Compiled defaults
      type(api_model_options_v1_0), target :: defaults
      !> Default values; only named components are copied to the wire layout
      type(api_model_options_v1_0) :: values
      !> Decoded error handle
      type(vp_error), pointer :: error
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call zero_options_storage(c_loc(defaults), c_sizeof(defaults))
      defaults%struct_size = bytes
      defaults%debug = values%debug
      defaults%verbosity = values%verbosity
      call copy_options(c_loc(defaults), options, c_sizeof(defaults), bytes, model_options_min_size, error%ptr)
      call prefix_api_error(error%ptr, "init_model_options")
   end subroutine init_model_options_api

   !> Decode optional model options into an independent value
   subroutine read_model_options(options, value, error)
      !> Optional caller options; NULL selects defaults
      type(c_ptr), intent(in) :: options
      !> Default-initialized result
      type(api_model_options_v1_0), intent(out), target :: value
      !> Diagnostic on invalid size
      type(error_type), allocatable, intent(out) :: error
      !> Caller allocation size, the common first field
      integer(c_size_t), pointer :: bytes
      value = api_model_options_v1_0()
      value%struct_size = c_sizeof(value)
      if (.not. c_associated(options)) return
      call c_f_pointer(options, bytes)
      call copy_options(options, c_loc(value), c_sizeof(value), bytes, model_options_min_size, error)
   end subroutine read_model_options

   !> Initialize pcm options within the caller's allocation
   subroutine init_pcm_options_api(verror, options, bytes) bind(C, name="moist_init_pcm_options")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Caller-owned options buffer
      type(c_ptr), value, intent(in) :: options
      !> Allocated byte count
      integer(c_size_t), value, intent(in) :: bytes
      !> Compiled defaults
      type(api_pcm_options_v1_0), target :: defaults
      !> Default values; only named components are copied to the wire layout
      type(api_pcm_options_v1_0) :: values
      !> Decoded error handle
      type(vp_error), pointer :: error
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      values = default_pcm_options()
      call zero_options_storage(c_loc(defaults), c_sizeof(defaults))
      defaults%struct_size = bytes
      defaults%solver = values%solver
      defaults%solver_maxiter = values%solver_maxiter
      defaults%solver_tol = values%solver_tol
      call copy_options(c_loc(defaults), options, c_sizeof(defaults), bytes, pcm_options_min_size, error%ptr)
      call prefix_api_error(error%ptr, "init_pcm_options")
   end subroutine init_pcm_options_api

   !> Decode optional pcm options into an independent value
   subroutine read_pcm_options(options, value, error)
      !> Optional caller options; NULL selects defaults
      type(c_ptr), intent(in) :: options
      !> Default-initialized result
      type(api_pcm_options_v1_0), intent(out), target :: value
      !> Diagnostic on invalid size
      type(error_type), allocatable, intent(out) :: error
      !> Caller allocation size, the common first field
      integer(c_size_t), pointer :: bytes
      value = default_pcm_options()
      value%struct_size = c_sizeof(value)
      if (.not. c_associated(options)) return
      call c_f_pointer(options, bytes)
      call copy_options(options, c_loc(value), c_sizeof(value), bytes, pcm_options_min_size, error)
   end subroutine read_pcm_options

   !> Clear all bytes before initializing individual fields, including ABI padding
   subroutine zero_options_storage(address, length)
      !> Options object to initialize
      type(c_ptr), intent(in) :: address
      !> Complete allocation size
      integer(c_size_t), intent(in) :: length
      !> Byte view, including padding between and after fields
      integer(c_int8_t), pointer :: bytes(:)

      call c_f_pointer(address, bytes, [length])
      bytes = 0_c_int8_t
   end subroutine zero_options_storage

   !> Copy the supported options prefix after validating its allocated size
   subroutine copy_options(source, destination, required, available, minimum, error)
      !> Source address
      type(c_ptr), intent(in) :: source
      !> Destination address
      type(c_ptr), intent(in) :: destination
      !> Supported prefix length
      integer(c_size_t), intent(in) :: required
      !> Caller allocation length
      integer(c_size_t), intent(in) :: available
      !> Frozen minimum size of the first supported ABI layout
      integer(c_size_t), intent(in) :: minimum
      !> Diagnostic; destination is unchanged on error
      type(error_type), allocatable, intent(out) :: error
      !> Byte views, bounded to the supported prefix
      integer(c_int8_t), pointer :: src(:), dst(:)
      !> Complete shared prefix length
      integer(c_size_t) :: copied
      if (.not. c_associated(source) .or. .not. c_associated(destination)) then
         call fatal_error(error, "Options pointer is missing")
         return
      end if
      if (available < minimum) then
         call fatal_error(error, "Options struct_size is smaller than the 1.0 layout")
         return
      end if
      copied = min(required, available)
      call c_f_pointer(source, src, [copied])
      call c_f_pointer(destination, dst, [copied])
      dst = src
   end subroutine copy_options

   !> Delete an independently owned LSF handle
   subroutine delete_lsf_api(handle) bind(C, name="moist_delete_lsf")
      !> Handle address, set to NULL on return
      type(c_ptr), intent(inout), optional :: handle
      !> Decoded handle
      type(vp_lsf), pointer :: lsf
      if (.not. present(handle)) return
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, lsf)
      deallocate (lsf)
      handle = c_null_ptr
   end subroutine delete_lsf_api

   !> Create an independent svdw level-set function
   function new_svdw_lsf_api(verror, options) result(handle) bind(C, name="moist_new_svdw_lsf")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Optional settings
      type(c_ptr), value, intent(in) :: options
      !> Owned LSF, NULL on failure
      type(c_ptr) :: handle
      type(vp_error), pointer :: error
      type(vp_lsf), pointer :: lsf
      type(api_svdw_options_v1_0) :: o
      integer :: stat
      nullify (lsf)
      handle = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call read_svdw_options(options, o, error%ptr)
      call prefix_api_error(error%ptr, "new_svdw_lsf")
      if (allocated(error%ptr)) return
      allocate (lsf, stat=stat)
      if (stat == 0) allocate (moist_cavity_drop_lsf_svdw_type :: lsf%ptr, stat=stat)
      if (stat /= 0) then
         if (associated(lsf)) deallocate (lsf)
         call api_error(error%ptr, "new_svdw_lsf", "Cannot allocate LSF")
         return
      end if
      select type (item => lsf%ptr)
      type is (moist_cavity_drop_lsf_svdw_type)
         call item%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=real(o%blend_k, wp), &
                                       blend_1b=real(o%blend_1b, wp), blend_2b=real(o%blend_2b, wp), blend_3b=real(o%blend_3b, wp)))
      end select
      handle = c_loc(lsf)
   end function new_svdw_lsf_api

   !> Create an independent cfc level-set function
   function new_cfc_lsf_api(verror, options) result(handle) bind(C, name="moist_new_cfc_lsf")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Optional settings
      type(c_ptr), value, intent(in) :: options
      !> Owned LSF, NULL on failure
      type(c_ptr) :: handle
      type(vp_error), pointer :: error
      type(vp_lsf), pointer :: lsf
      type(api_cfc_options_v1_0) :: o
      integer :: stat
      nullify (lsf)
      handle = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call read_cfc_options(options, o, error%ptr)
      call prefix_api_error(error%ptr, "new_cfc_lsf")
      if (allocated(error%ptr)) return
      allocate (lsf, stat=stat)
      if (stat == 0) allocate (moist_cavity_drop_lsf_cfc_type :: lsf%ptr, stat=stat)
      if (stat /= 0) then
         if (associated(lsf)) deallocate (lsf)
         call api_error(error%ptr, "new_cfc_lsf", "Cannot allocate LSF")
         return
      end if
      select type (item => lsf%ptr)
      type is (moist_cavity_drop_lsf_cfc_type)
         call item%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=real(o%a1, wp), a2=real(o%a2, wp), &
                                                                  c=real(o%c, wp), m=int(o%m)))
      end select
      handle = c_loc(lsf)
   end function new_cfc_lsf_api

   !> Create an LSF borrowing a host callback and its context
   function new_isodensity_callback_lsf_api(verror, callback, context, options) result(handle) &
      bind(C, name="moist_new_isodensity_callback_lsf")
      !> Diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Host density callback, borrowed for the lifetime of the LSF
      type(c_funptr), value, intent(in) :: callback
      !> Opaque context passed to the host callback
      type(c_ptr), value, intent(in) :: context
      !> Optional settings; NULL selects compiled defaults
      type(c_ptr), value, intent(in) :: options
      type(c_ptr) :: handle
      type(vp_error), pointer :: error
      type(vp_lsf), pointer :: lsf
      type(api_isodensity_options_v1_0) :: o
      integer :: stat
      nullify (lsf)
      handle = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call read_isodensity_options(options, o, error%ptr)
      call prefix_api_error(error%ptr, "new_isodensity_callback_lsf")
      if (allocated(error%ptr)) return
      if (.not. c_associated(callback)) then
         call api_error(error%ptr, "new_isodensity_callback_lsf", "Callback is missing")
         return
      end if
      if (.not. valid_isodensity_options(o, error%ptr)) then
         call prefix_api_error(error%ptr, "new_isodensity_callback_lsf")
         return
      end if
      allocate (lsf, stat=stat)
      if (stat == 0) allocate (moist_cavity_drop_lsf_isodensity_callback_type :: lsf%ptr, stat=stat)
      if (stat /= 0) then
         if (associated(lsf)) deallocate (lsf)
         call api_error(error%ptr, "new_isodensity_callback_lsf", "Cannot allocate LSF")
         return
      end if
      select type (item => lsf%ptr)
      type is (moist_cavity_drop_lsf_isodensity_callback_type)
         call item%new(callback_ptr=callback, context=context, &
                       param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=real(o%rho_iso, wp), &
                                                                         scale=real(o%scale, wp)))
      end select
      handle = c_loc(lsf)
   end function new_isodensity_callback_lsf_api

   !> Validate the physical isodensity controls
   function valid_isodensity_options(o, error) result(valid)
      !> Isodensity settings to validate
      type(api_isodensity_options_v1_0), intent(in) :: o
      !> Diagnostic on failure
      type(error_type), allocatable, intent(out) :: error
      logical :: valid
      valid = o%rho_iso > 0.0_c_double .and. o%scale > 0.0_c_double
      if (.not. valid) call fatal_error(error, "rho_iso and scale must be positive")
   end function valid_isodensity_options

   !> Create an LSF with its own Cartesian Gaussian basis
   function new_isodensity_lsf_api(verror, nshell, c_shell_atom, c_shell_l, &
                                   c_shell_nprim, c_exps, c_coeffs, options) result(handle) bind(C, name="moist_new_isodensity_lsf")
      !> Diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Number of contracted Gaussian shells
      integer(c_int), value, intent(in) :: nshell
      !> Zero-based atom index for each shell
      type(c_ptr), value, intent(in) :: c_shell_atom
      !> Angular momentum of each shell
      type(c_ptr), value, intent(in) :: c_shell_l
      !> Primitive count for each shell
      type(c_ptr), value, intent(in) :: c_shell_nprim
      !> Packed primitive exponents, grouped by shell
      type(c_ptr), value, intent(in) :: c_exps
      !> Packed contraction coefficients, including primitive normalization
      type(c_ptr), value, intent(in) :: c_coeffs
      !> Optional settings; NULL selects compiled defaults
      type(c_ptr), value, intent(in) :: options
      type(c_ptr) :: handle
      type(vp_error), pointer :: error
      type(vp_lsf), pointer :: lsf
      type(api_isodensity_options_v1_0) :: o
      integer(c_int), pointer :: atoms(:), angular(:), primitives(:)
      real(c_double), pointer :: exps(:), coeffs(:)
      integer :: nprim, stat
      nullify (lsf)
      handle = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call read_isodensity_options(options, o, error%ptr)
      call prefix_api_error(error%ptr, "new_isodensity_lsf")
      if (allocated(error%ptr)) return
      if (.not. valid_isodensity_options(o, error%ptr)) then
         call prefix_api_error(error%ptr, "new_isodensity_lsf")
         return
      end if
      if (nshell < 1 .or. .not. c_associated(c_shell_atom) .or. .not. c_associated(c_shell_l) &
          .or. .not. c_associated(c_shell_nprim) .or. .not. c_associated(c_exps) &
          .or. .not. c_associated(c_coeffs)) then
         call api_error(error%ptr, "new_isodensity_lsf", "Missing basis arrays or invalid shell count")
         return
      end if
      call c_f_pointer(c_shell_atom, atoms, [nshell])
      call c_f_pointer(c_shell_l, angular, [nshell])
      call c_f_pointer(c_shell_nprim, primitives, [nshell])
      if (any(primitives < 1) .or. sum(int(primitives, c_int64_t)) > huge(nprim)) then
         call api_error(error%ptr, "new_isodensity_lsf", "Invalid primitive counts")
         return
      end if
      nprim = sum(primitives)
      call c_f_pointer(c_exps, exps, [nprim])
      call c_f_pointer(c_coeffs, coeffs, [nprim])
      allocate (lsf, stat=stat)
      if (stat == 0) allocate (moist_cavity_drop_lsf_isodensity_internal_type :: lsf%ptr, stat=stat)
      if (stat /= 0) then
         if (associated(lsf)) deallocate (lsf)
         call api_error(error%ptr, "new_isodensity_lsf", "Cannot allocate LSF")
         return
      end if
      select type (item => lsf%ptr)
      type is (moist_cavity_drop_lsf_isodensity_internal_type)
         call item%new(sh_atom=int(atoms) + 1, sh_l=int(angular), sh_nprim=int(primitives), exps=real(exps, wp), &
                       coeffs=real(coeffs, wp), error=error%ptr, &
                       param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=real(o%rho_iso, wp), &
                                                                         scale=real(o%scale, wp)))
      end select
      call prefix_api_error(error%ptr, "new_isodensity_lsf")
      if (allocated(error%ptr)) then
         deallocate (lsf)
         return
      end if
      handle = c_loc(lsf)
   end function new_isodensity_lsf_api

   !> Resolve optional radii; NULL selects CPCM radii
   subroutine read_api_radii(handle, radii, error)
      !> Optional radii handle; NULL selects CPCM radii
      type(c_ptr), intent(in) :: handle
      !> Independent copy of the selected radii model
      class(radius_type), allocatable, intent(out) :: radii
      !> Diagnostic on failure
      type(error_type), allocatable, intent(out) :: error
      type(vp_radii), pointer :: source
      integer :: stat
      if (.not. c_associated(handle)) then
         call new_radii("cpcm", radii, error)
         return
      end if
      call c_f_pointer(handle, source)
      if (.not. allocated(source%ptr)) then
         call fatal_error(error, "Radii are not initialized")
         return
      end if
      allocate (radii, source=source%ptr, stat=stat)
      if (stat /= 0) call fatal_error(error, "Cannot copy radii")
   end subroutine read_api_radii

   !> Construct a drop cavity from copied configuration
   function create_drop_cavity_api(verror, vlsf, vradii, options) result(handle) bind(C, name="moist_new_drop_cavity")
      !> Diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Level-set function handle, copied into the cavity
      type(c_ptr), value, intent(in) :: vlsf
      type(vp_lsf), pointer :: lsf
      !> Radii handle; NULL selects defaults in constructors
      type(c_ptr), value, intent(in) :: vradii
      !> Optional settings; NULL selects compiled defaults
      type(c_ptr), value, intent(in) :: options
      type(c_ptr) :: handle
      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav
      type(api_drop_options_v1_0) :: o
      class(radius_type), allocatable :: radii
      integer :: stat
      nullify (lsf, cav)
      handle = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call read_drop_options(options, o, error%ptr)
      call prefix_api_error(error%ptr, "new_drop_cavity")
      if (allocated(error%ptr)) return
      if (.not. c_associated(vlsf)) then
         call api_error(error%ptr, "new_drop_cavity", "LSF handle is required")
         return
      end if
      call c_f_pointer(vlsf, lsf)
      if (.not. allocated(lsf%ptr)) then
         call api_error(error%ptr, "new_drop_cavity", "LSF is not initialized")
         return
      end if
      if (o%tolerance <= 0.0_c_double .or. o%rho_grid_h <= 0.0_c_double .or. &
          o%branch_weight_s < 0.0_c_double .or. o%proj_maxiter < 1 .or. o%proj_level < 1 .or. o%proj_level > 9) then
         call api_error(error%ptr, "new_drop_cavity", "Invalid DROP options")
         return
      end if
      call read_api_radii(vradii, radii, error%ptr)
      call prefix_api_error(error%ptr, "new_drop_cavity")
      if (allocated(error%ptr)) return
      allocate (cav, stat=stat)
      if (stat == 0) allocate (cavity_type_drop :: cav%ptr, stat=stat)
      if (stat /= 0) then
         if (associated(cav)) deallocate (cav)
         call api_error(error%ptr, "new_drop_cavity", "Cannot allocate cavity")
         return
      end if
      call new_context(cav%ctx, verbosity=int(o%verbosity), debug=logical(o%debug))
      select type (item => cav%ptr)
      type is (cavity_type_drop)
         call new_cavity_drop(item, cav%ctx, radius_model=radii, error=error%ptr, lsf_model=lsf%ptr, &
                              param=moist_cavity_drop_parameters_type(num_leb=int(o%nleb), tolerance=real(o%tolerance, wp), &
                                       do_fine=logical(o%do_fine), proj_maxiter=int(o%proj_maxiter), proj_level=int(o%proj_level), &
                                                   branch_weight_s=real(o%branch_weight_s, wp), rho_grid_h=real(o%rho_grid_h, wp), &
                                                                      wleb_prune_level=int(o%wleb_prune_level)))
      end select
      call prefix_api_error(error%ptr, "new_drop_cavity")
      if (allocated(error%ptr)) then
         if (associated(cav%ptr)) deallocate (cav%ptr)
         call cav%ctx%delete()
         deallocate (cav)
         return
      end if
      handle = c_loc(cav)
   end function create_drop_cavity_api

   !> Construct a iswig cavity from copied configuration
   function create_iswig_cavity_api(verror, vradii, options) result(handle) bind(C, name="moist_new_iswig_cavity")
      !> Diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Radii handle; NULL selects defaults in constructors
      type(c_ptr), value, intent(in) :: vradii
      !> Optional settings; NULL selects compiled defaults
      type(c_ptr), value, intent(in) :: options
      type(c_ptr) :: handle
      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav
      type(api_iswig_options_v1_0) :: o
      class(radius_type), allocatable :: radii
      integer :: stat
      nullify (cav)
      handle = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call read_iswig_options(options, o, error%ptr)
      call prefix_api_error(error%ptr, "new_iswig_cavity")
      if (allocated(error%ptr)) return
      call read_api_radii(vradii, radii, error%ptr)
      call prefix_api_error(error%ptr, "new_iswig_cavity")
      if (allocated(error%ptr)) return
      allocate (cav, stat=stat)
      if (stat == 0) allocate (cavity_type_iswig :: cav%ptr, stat=stat)
      if (stat /= 0) then
         if (associated(cav)) deallocate (cav)
         call api_error(error%ptr, "new_iswig_cavity", "Cannot allocate cavity")
         return
      end if
      call new_context(cav%ctx, verbosity=int(o%verbosity), debug=logical(o%debug))
      select type (item => cav%ptr)
      type is (cavity_type_iswig)
         call new_cavity_iswig(item, cav%ctx, radius_model=radii, error=error%ptr, &
                               param=moist_cavity_iswig_parameters_type(num_leb=int(o%nleb), cut_a=real(o%cut_a, wp), &
                                                                        cut_f=real(o%cut_f, wp)))
      end select
      call prefix_api_error(error%ptr, "new_iswig_cavity")
      if (allocated(error%ptr)) then
         if (associated(cav%ptr)) deallocate (cav%ptr)
         call cav%ctx%delete()
         deallocate (cav)
         return
      end if
      handle = c_loc(cav)
   end function create_iswig_cavity_api

   !> Create a model with an owned cavity copy and optional logging settings
   function create_model_api(verror, cavity, options) result(handle) bind(C, name="moist_new_model")
      !> Diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Cavity handle, copied into the model
      type(c_ptr), value, intent(in) :: cavity
      !> Optional settings; NULL selects compiled defaults
      type(c_ptr), value, intent(in) :: options
      type(c_ptr) :: handle
      type(vp_error), pointer :: error
      type(api_model_options_v1_0) :: o
      handle = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call read_model_options(options, o, error%ptr)
      call prefix_api_error(error%ptr, "new_model")
      if (allocated(error%ptr)) return
      handle = new_general_solvation_model_api(verror, cavity, o%debug, o%verbosity)
   end function create_model_api

   !> Create a cpcm component with optional solver settings
   function create_cpcm_component_api(verror, epsilon, options) result(handle) bind(C, name="moist_new_cpcm_component")
      !> Diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Relative dielectric constant
      real(c_double), value, intent(in) :: epsilon
      !> Optional settings; NULL selects compiled defaults
      type(c_ptr), value, intent(in) :: options
      type(c_ptr) :: handle
      type(vp_error), pointer :: error
      type(api_pcm_options_v1_0) :: o
      handle = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call read_pcm_options(options, o, error%ptr)
      call prefix_api_error(error%ptr, "new_cpcm_component")
      if (allocated(error%ptr)) return
      handle = new_cpcm_component_api(verror, epsilon, o%solver, o%solver_tol, o%solver_maxiter)
   end function create_cpcm_component_api

   !> Create a cosmo component with optional solver settings
   function create_cosmo_component_api(verror, epsilon, options) result(handle) bind(C, name="moist_new_cosmo_component")
      !> Diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Relative dielectric constant
      real(c_double), value, intent(in) :: epsilon
      !> Optional settings; NULL selects compiled defaults
      type(c_ptr), value, intent(in) :: options
      type(c_ptr) :: handle
      type(vp_error), pointer :: error
      type(api_pcm_options_v1_0) :: o
      handle = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call read_pcm_options(options, o, error%ptr)
      call prefix_api_error(error%ptr, "new_cosmo_component")
      if (allocated(error%ptr)) return
      handle = new_cosmo_component_api(verror, epsilon, o%solver, o%solver_tol, o%solver_maxiter)
   end function create_cosmo_component_api

   !> Build a diagnostic tagged with its calling routine
   subroutine api_error(error, routine, msg)
      !> Diagnostic on failure
      type(error_type), allocatable, intent(out) :: error
      !> Public C entry point, without the common moist_ prefix
      character(len=*), intent(in) :: routine
      !> Diagnostic details
      character(len=*), intent(in) :: msg
      call fatal_error(error, "["//namespace//routine//"] "//msg)
   end subroutine api_error

   !> Add the calling routine to a propagated error
   subroutine prefix_api_error(error, routine)
      !> Error to annotate
      type(error_type), allocatable, intent(inout) :: error
      !> Calling entry point
      character(len=*), intent(in) :: routine
      !> Saved diagnostic before replacing the error
      character(len=:), allocatable :: message
      if (.not. allocated(error)) return
      message = error%message
      ! Replace an inner C entry point's prefix when a public routine delegates
      if (index(message, "["//namespace) == 1) then
         message = message(index(message, "] ") + 2:)
      end if
      call api_error(error, routine, message)
   end subroutine prefix_api_error

   !> Obtain library version as major * 10000 + minor * 100 + patch
   function get_version_api() result(version) &
         & bind(C, name=namespace//"get_version")
      integer(c_int) :: version
      integer :: major, minor, patch

      call get_moist_version(major, minor, patch)
      version = 10000_c_int*major + 100_c_int*minor + patch

   end function get_version_api

   !> Copy the full release version or query its length, including the prerelease suffix
   subroutine get_version_string_api(verror, buffer, capacity, length) bind(C, name="moist_get_version_string")
      !> Required diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Output buffer; NULL is allowed only with zero capacity
      character(kind=c_char), intent(inout), optional :: buffer(*)
      !> Buffer capacity in bytes, including the terminator
      integer(c_size_t), value, intent(in) :: capacity
      !> Required text length excluding the terminator; unchanged on error
      integer(c_size_t), intent(inout), optional :: length
      !> Full release version
      character(len=:), allocatable :: version
      !> Decoded diagnostic handle
      type(vp_error), pointer :: error

      if (.not. valid_string_output(verror, buffer, capacity, length, "get_version_string", error)) return
      call get_moist_version(string=version)
      call copy_string_output(version, buffer, capacity, length)
   end subroutine get_version_string_api

   !> Create an error handle
   function new_error_api() &
         & result(verror) &
         & bind(C, name=namespace//"new_error")
      type(vp_error), pointer :: error
      type(c_ptr) :: verror
      !> Allocation status; no diagnostic handle exists yet
      integer :: stat

      verror = c_null_ptr
      allocate (error, stat=stat)
      if (stat == 0) verror = c_loc(error)

   end function new_error_api

   !> Delete an error handle
   subroutine delete_error_api(verror) &
         & bind(C, name=namespace//"delete_error")
      !> Diagnostic handle
      type(c_ptr), intent(inout), optional :: verror
      type(vp_error), pointer :: error

      if (.not. present(verror)) return
      if (c_associated(verror)) then
         call c_f_pointer(verror, error)

         deallocate (error)
         verror = c_null_ptr
      end if

   end subroutine delete_error_api

   !> Query error handle status
   function check_error_api(verror) result(status) &
         & bind(C, name=namespace//"check_error")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      integer(c_int) :: status

      if (c_associated(verror)) then
         call c_f_pointer(verror, error)

         if (allocated(error%ptr)) then
            status = api_failure
         else
            status = api_success
         end if
      else
         status = api_invalid_error
      end if

   end function check_error_api

   !> Copy a bounded diagnostic without modifying the error handle
   subroutine get_error_api(verror, charptr, buffersize) bind(C, name="moist_get_error")
      !> Diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Diagnostic buffer; always null-terminated when capacity permits
      character(kind=c_char), intent(inout), optional :: charptr(*)
      !> Allocated buffer size in bytes
      integer(c_int), intent(in), optional :: buffersize
      type(vp_error), pointer :: error
      if (.not. present(charptr) .or. .not. present(buffersize)) return
      if (buffersize < 1) return
      charptr(1) = c_null_char
      if (.not. c_associated(verror)) then
         call f_c_character("[moist_get_error] Invalid error handle", charptr, buffersize)
         return
      end if
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) call f_c_character(error%ptr%message, charptr, buffersize)
   end subroutine get_error_api

   !> Create new molecular structure data (quantities in Bohr)
   function new_structure_api(verror, natoms, numbers, positions, &
         & c_lattice, c_periodic) result(vmol) &
         & bind(C, name=namespace//"new_structure")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Number of atoms
      integer(c_int), value, intent(in) :: natoms
      !> Atomic numbers
      integer(c_int), intent(in), optional :: numbers(natoms)
      !> Atomic positions in Bohr, shape (3,natoms)
      real(c_double), intent(in), optional :: positions(3, natoms)
      !> Optional lattice vectors in Bohr
      real(c_double), intent(in), optional :: c_lattice(3, 3)
      !> Optional periodicity flags by lattice direction
      logical(c_bool), intent(in), optional :: c_periodic(3)
      !> Owned result, NULL on failure
      type(c_ptr) :: vmol
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      real(wp), allocatable :: lattice(:, :)
      logical, allocatable :: periodic(:)
      type(vp_structure), pointer :: mol

      vmol = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(numbers)) then
         call api_error(error%ptr, "new_structure", "Required pointer 'numbers' is missing")
         return
      end if
      if (.not. present(positions)) then
         call api_error(error%ptr, "new_structure", "Required pointer 'positions' is missing")
         return
      end if
      if (present(c_lattice)) then
         allocate (lattice(3, 3))
         lattice(:, :) = c_lattice
      end if
      if (present(c_periodic)) then
         allocate (periodic(3))
         periodic(:) = c_periodic
      end if

      allocate (mol)
      call new(mol%ptr, numbers, positions, lattice=lattice, periodic=periodic)
      vmol = c_loc(mol)

      call verify_structure(error%ptr, mol%ptr)
      call prefix_api_error(error%ptr, "new_structure")

   end function new_structure_api

   !> Delete molecular structure data
   subroutine delete_structure_api(vmol) &
         & bind(C, name=namespace//"delete_structure")
      !> Molecular structure handle
      type(c_ptr), intent(inout), optional :: vmol
      type(vp_structure), pointer :: mol

      if (.not. present(vmol)) return
      if (c_associated(vmol)) then
         call c_f_pointer(vmol, mol)

         deallocate (mol)
         vmol = c_null_ptr
      end if

   end subroutine delete_structure_api

   !> Update coordinates and lattice parameters (quantities in Bohr)
   subroutine update_structure_api(verror, vmol, positions, lattice) &
         & bind(C, name=namespace//"update_structure")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Molecular structure handle
      type(c_ptr), value :: vmol
      !> Atomic positions in Bohr, shape (3,natoms)
      real(c_double), intent(in), optional :: positions(3, *)
      !> Lattice vectors in Bohr
      real(c_double), intent(in), optional :: lattice(3, 3)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      type(vp_structure), pointer :: mol

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(positions)) then
         call api_error(error%ptr, "update_structure", "Required pointer 'positions' is missing")
         return
      end if
      if (.not. c_associated(vmol)) then
         call api_error(error%ptr, "update_structure", "Molecular structure data is missing")
         return
      end if
      call c_f_pointer(vmol, mol)

      if (mol%ptr%nat <= 0 .or. mol%ptr%nid <= 0 .or. .not. allocated(mol%ptr%num) &
         & .or. .not. allocated(mol%ptr%id) .or. .not. allocated(mol%ptr%xyz)) then
         call api_error(error%ptr, "update_structure", "Invalid molecular structure data provided")
         return
      end if

      mol%ptr%xyz(:, :) = positions(:3, :mol%ptr%nat)
      if (present(lattice)) then
         mol%ptr%lattice(:, :) = lattice(:3, :3)
      end if

      call verify_structure(error%ptr, mol%ptr)
      call prefix_api_error(error%ptr, "update_structure")

   end subroutine update_structure_api

   !> Create a new radii handle from a named model
   subroutine new_radii_handle_api(verror, model_name, routine_name, vradii)
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Built-in radii model name
      character(len=*), intent(in) :: model_name
      !> Public C entry point, without the common moist_ prefix
      character(len=*), intent(in) :: routine_name
      !> Radii handle, created or released by this call
      type(c_ptr), intent(out) :: vradii
      type(vp_radii), pointer :: radii
      type(error_type), allocatable :: radii_error

      vradii = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      allocate (radii)
      call new_radii(model_name, radii%ptr, radii_error)
      if (allocated(radii_error)) then
         call api_error(error%ptr, routine_name, radii_error%message)
         if (allocated(radii%ptr)) deallocate (radii%ptr)
         deallocate (radii)
         return
      end if

      vradii = c_loc(radii)
   end subroutine new_radii_handle_api

   !> Create CPCM radii model handle
   function new_cpcm_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_cpcm_radii")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "cpcm", "new_cpcm_radii", vradii)
   end function new_cpcm_radii_api

   !> Create SMD radii model handle
   function new_smd_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_smd_radii")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "smd", "new_smd_radii", vradii)
   end function new_smd_radii_api

   !> Create D3 radii model handle
   function new_d3_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_d3_radii")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "d3", "new_d3_radii", vradii)
   end function new_d3_radii_api

   !> Create COSMO radii model handle
   function new_cosmo_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_cosmo_radii")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "cosmo", "new_cosmo_radii", vradii)
   end function new_cosmo_radii_api

   !> Create Bondi radii model handle
   function new_bondi_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_bondi_radii")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "bondi", "new_bondi_radii", vradii)
   end function new_bondi_radii_api

   !> Create custom radii model handle
   !>
   !> - initialize with set_custom_radii_atoms or set_custom_radii_elements before use
   function new_custom_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_custom_radii")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(vp_radii), pointer :: radii
      type(c_ptr) :: vradii

      vradii = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      allocate (radii)
      allocate (radius_type_custom :: radii%ptr)
      vradii = c_loc(radii)
   end function new_custom_radii_api

   !> Set custom radii from per-atom values (bohr)
   subroutine set_custom_radii_atoms_api(verror, vradii, natoms, atom_radii) &
         & bind(C, name=namespace//"set_custom_radii_atoms")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Radii model handle
      type(c_ptr), value :: vradii
      !> Number of atoms
      integer(c_int), value :: natoms
      !> Per-atom radii in Bohr
      real(c_double), intent(in), optional :: atom_radii(*)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      type(vp_radii), pointer :: radii
      type(error_type), allocatable :: radii_error
      real(wp), allocatable :: atom_radii_wp(:)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(atom_radii)) then
         call api_error(error%ptr, "set_custom_radii_atoms", "Required pointer 'atom_radii' is missing")
         return
      end if
      if (.not. c_associated(vradii)) then
         call api_error(error%ptr, "set_custom_radii_atoms", "Radii handle is missing")
         return
      end if
      call c_f_pointer(vradii, radii)

      if (.not. allocated(radii%ptr)) then
         call api_error(error%ptr, "set_custom_radii_atoms", "Radii model is not initialized")
         return
      end if

      if (natoms < 1) then
         call api_error(error%ptr, "set_custom_radii_atoms", "natoms must be positive")
         return
      end if

      allocate (atom_radii_wp(natoms))
      atom_radii_wp(:) = atom_radii(:natoms)

      select type (model => radii%ptr)
      type is (radius_type_custom)
         call new_custom_radii_atoms(atom_radii_wp, model, radii_error)
         if (allocated(radii_error)) then
            call api_error(error%ptr, "set_custom_radii_atoms", radii_error%message)
         end if
      class default
         call api_error(error%ptr, "set_custom_radii_atoms", "Radii model is not custom type")
      end select

   end subroutine set_custom_radii_atoms_api

   !> Set custom radii from per-element values (bohr)
   subroutine set_custom_radii_elements_api(verror, vradii, nentries, atomic_numbers, element_radii) &
         & bind(C, name=namespace//"set_custom_radii_elements")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Radii model handle
      type(c_ptr), value :: vradii
      !> Number of supplied element radii
      integer(c_int), value :: nentries
      !> Atomic numbers
      integer(c_int), intent(in), optional :: atomic_numbers(*)
      !> Per-element radii in Bohr
      real(c_double), intent(in), optional :: element_radii(*)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      type(vp_radii), pointer :: radii
      type(error_type), allocatable :: radii_error
      integer, allocatable :: atomic_numbers_f(:)
      real(wp), allocatable :: element_radii_wp(:)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(atomic_numbers)) then
         call api_error(error%ptr, "set_custom_radii_elements", "Required pointer 'atomic_numbers' is missing")
         return
      end if
      if (.not. present(element_radii)) then
         call api_error(error%ptr, "set_custom_radii_elements", "Required pointer 'element_radii' is missing")
         return
      end if
      if (.not. c_associated(vradii)) then
         call api_error(error%ptr, "set_custom_radii_elements", "Radii handle is missing")
         return
      end if
      call c_f_pointer(vradii, radii)

      if (.not. allocated(radii%ptr)) then
         call api_error(error%ptr, "set_custom_radii_elements", "Radii model is not initialized")
         return
      end if

      if (nentries < 1) then
         call api_error(error%ptr, "set_custom_radii_elements", "nentries must be positive")
         return
      end if

      allocate (atomic_numbers_f(nentries), element_radii_wp(nentries))
      atomic_numbers_f(:) = atomic_numbers(:nentries)
      element_radii_wp(:) = element_radii(:nentries)

      select type (model => radii%ptr)
      type is (radius_type_custom)
         call new_custom_radii_elements(atomic_numbers_f, element_radii_wp, model, radii_error)
         if (allocated(radii_error)) then
            call api_error(error%ptr, "set_custom_radii_elements", radii_error%message)
         end if
      class default
         call api_error(error%ptr, "set_custom_radii_elements", "Radii model is not custom type")
      end select

   end subroutine set_custom_radii_elements_api

   !> Delete radii model handle
   subroutine delete_radii_api(vradii) &
         & bind(C, name=namespace//"delete_radii")
      !> Radii handle, created or released by this call
      type(c_ptr), intent(inout), optional :: vradii
      type(vp_radii), pointer :: radii

      if (.not. present(vradii)) return
      if (c_associated(vradii)) then
         call c_f_pointer(vradii, radii)
         if (allocated(radii%ptr)) deallocate (radii%ptr)
         deallocate (radii)
         vradii = c_null_ptr
      end if
   end subroutine delete_radii_api

   !> Delete solvation model handle
   subroutine delete_solvation_model_api(vmodel) &
         & bind(C, name=namespace//"delete_model")
      !> Model handle
      type(c_ptr), intent(inout), optional :: vmodel
      type(vp_model), pointer :: model

      if (.not. present(vmodel)) return
      if (c_associated(vmodel)) then
         call c_f_pointer(vmodel, model)

         if (allocated(model%ptr)) deallocate (model%ptr)
         call model%ctx%delete()
         deallocate (model)
         vmodel = c_null_ptr
      end if

   end subroutine delete_solvation_model_api

   !> Update a solvation model with a molecular structure
   subroutine update_solvation_model_api(verror, vmodel, vmol) &
         & bind(C, name=namespace//"update_model")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Model handle
      type(c_ptr), value :: vmodel
      type(vp_model), pointer :: model
      !> Molecular structure handle
      type(c_ptr), value :: vmol
      type(vp_structure), pointer :: mol
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vmodel)) then
         call api_error(error%ptr, "update_model", "Model handle is missing")
         return
      end if
      call c_f_pointer(vmodel, model)

      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "update_model", "Model is not initialized")
         return
      end if

      if (.not. c_associated(vmol)) then
         call api_error(error%ptr, "update_model", "Molecular structure data is missing")
         return
      end if
      call c_f_pointer(vmol, mol)

      call model%ptr%update(mol%ptr, error=model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "update_model", model_error%message)
         return
      end if

   end subroutine update_solvation_model_api

   !> Get a borrowed cavity handle from a solvation model
   !>
   !> - NOT owned by the caller; independent cavity-update entry points reject it
   !> - moist_delete_cavity releases only the borrowed wrapper
   !> - valid as long as the parent model exists
   function get_solvation_model_cavity_api(verror, vmodel) result(vcav) &
         & bind(C, name=namespace//"get_model_cavity")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Model handle
      type(c_ptr), value :: vmodel
      type(vp_model), pointer :: model
      type(c_ptr) :: vcav
      type(vp_cavity), pointer :: cav
      class(cavity_type), pointer :: cavity_ptr
      character(len=:), allocatable :: message

      vcav = c_null_ptr

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vmodel)) then
         call api_error(error%ptr, "get_model_cavity", "Model handle is missing")
         return
      end if
      call c_f_pointer(vmodel, model)

      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "get_model_cavity", "Model is not initialized")
         return
      end if

      call borrow_general_cavity(model%ptr, cavity_ptr, message)
      if (allocated(message)) then
         call api_error(error%ptr, "get_model_cavity", message)
         return
      end if

      allocate (cav)
      cav%ptr => cavity_ptr
      cav%owned = .false.
      vcav = c_loc(cav)

   end function get_solvation_model_cavity_api

   !> Point at the cavity a general solvation model owns, without taking it
   subroutine borrow_general_cavity(model, cavity_ptr, message)
      !> Solvation model that may own a cavity; borrowed pointer stays live,
      !> feeding a read-write handle
      class(solvation_model_type), intent(inout), target :: model
      !> Borrowed cavity, left null unless the model exposes one
      class(cavity_type), pointer, intent(out) :: cavity_ptr
      !> Reason the cavity could not be borrowed, unallocated on success
      character(len=:), allocatable, intent(out) :: message

      cavity_ptr => null()

      select type (general => model)
      type is (solvation_model_general)
         if (.not. allocated(general%cavity)) then
            message = "General model cavity is not initialized"
            return
         end if
         cavity_ptr => general%cavity
      class default
         message = "This solvation model type does not expose a cavity"
      end select

   end subroutine borrow_general_cavity

   !> Allocate either PCM-family component behind the common opaque handle
   subroutine new_pcm_component_common(verror, epsilon, solver, solver_tol, solver_maxiter, &
                                       use_cosmo, routine_name, vcomponent)
      !> Error handle
      type(c_ptr), value :: verror
      !> Relative dielectric constant
      real(c_double), value :: epsilon
      !> PCM solver enumeration
      integer(c_int), value :: solver
      !> Iterative solver residual threshold
      real(c_double), value :: solver_tol
      !> Iterative solver iteration cap
      integer(c_int), value :: solver_maxiter
      !> Select COSMO rather than CPCM
      logical, intent(in) :: use_cosmo
      !> Public routine name used in error messages
      character(len=*), intent(in) :: routine_name
      !> New component handle
      type(c_ptr), intent(out) :: vcomponent
      !> Decoded error handle
      type(vp_error), pointer :: error
      !> Component wrapper
      type(vp_component), pointer :: component
      !> Concrete PCM-family component
      class(solvation_model_component_type), allocatable :: item
      !> Constructor error
      type(error_type), allocatable :: component_error

      vcomponent = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      allocate (component)
      call new_context(component%ctx, verbosity=0, debug=.false.)
      if (use_cosmo) then
         allocate (solvation_model_component_cosmo :: item)
         select type (pcm => item)
         type is (solvation_model_component_cosmo)
            call new_component_cosmo(pcm, component%ctx, epsilon=real(epsilon, wp), error=component_error, &
                                     param=moist_pcm_parameters_type(solver=int(solver), &
                                     & solver_tol=real(solver_tol, wp), solver_maxiter=int(solver_maxiter)))
         end select
      else
         allocate (solvation_model_component_cpcm :: item)
         select type (pcm => item)
         type is (solvation_model_component_cpcm)
            call new_component_cpcm(pcm, component%ctx, epsilon=real(epsilon, wp), error=component_error, &
                                    param=moist_pcm_parameters_type(solver=int(solver), &
                                    & solver_tol=real(solver_tol, wp), solver_maxiter=int(solver_maxiter)))
         end select
      end if
      if (allocated(component_error)) then
         call api_error(error%ptr, routine_name, component_error%message)
         call component%ctx%delete()
         deallocate (component)
         return
      end if
      allocate (component%ptr, source=item)
      vcomponent = c_loc(component)

   end subroutine new_pcm_component_common

   !> Create a CPCM component handle for use with a general model
   function new_cpcm_component_api(verror, epsilon, solver, solver_tol, solver_maxiter) result(vcomponent)
      !> Diagnostic handle
      type(c_ptr), value :: verror
      !> Relative dielectric constant
      real(c_double), value :: epsilon
      !> PCM linear-solver selector
      integer(c_int), value :: solver
      !> Iterative solver residual threshold
      real(c_double), value :: solver_tol
      !> Iterative solver iteration cap
      integer(c_int), value :: solver_maxiter
      type(c_ptr) :: vcomponent

      call new_pcm_component_common(verror, epsilon, solver, solver_tol, solver_maxiter, .false., &
                                    "new_cpcm_component", vcomponent)

   end function new_cpcm_component_api

   !> Create a COSMO component handle for use with a general model
   function new_cosmo_component_api(verror, epsilon, solver, solver_tol, solver_maxiter) result(vcomponent)
      !> Diagnostic handle
      type(c_ptr), value :: verror
      !> Relative dielectric constant
      real(c_double), value :: epsilon
      !> PCM linear-solver selector
      integer(c_int), value :: solver
      !> Iterative solver residual threshold
      real(c_double), value :: solver_tol
      !> Iterative solver iteration cap
      integer(c_int), value :: solver_maxiter
      type(c_ptr) :: vcomponent

      call new_pcm_component_common(verror, epsilon, solver, solver_tol, solver_maxiter, .true., &
                                    "new_cosmo_component", vcomponent)

   end function new_cosmo_component_api

   !> Create a pressure-volume energy component handle
   function new_pv_component_api(verror, pressure) result(vcomponent) &
         & bind(C, name=namespace//"new_pv_component")
      !> Error handle
      type(c_ptr), value :: verror
      !> Pressure multiplying the cavity volume
      real(c_double), value :: pressure
      !> New component handle
      type(c_ptr) :: vcomponent
      !> Decoded error handle
      type(vp_error), pointer :: error
      !> Component wrapper
      type(vp_component), pointer :: component
      !> Concrete pressure-volume component
      type(solvation_model_component_pv) :: item

      vcomponent = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      allocate (component)
      call new_context(component%ctx, verbosity=0, debug=.false.)
      call new_component_pv(item, real(pressure, wp))
      item%ctx => component%ctx
      allocate (component%ptr, source=item)
      vcomponent = c_loc(component)

   end function new_pv_component_api

   !> Create a GOSTSHYP hydrostatic-pressure component handle
   !>
   !> - no own density traces; the host answers the coupling's Gaussian-moment
   !>   request (`moist_answer_coupling_request` with "gt", "pt", "mt", "rt") in every phase
   !> - amplitudes read back as the "gostshyp_amplitude" item of the response
   !>   walk (`moist_next_response_item`, then `moist_get_response_array`)
   function new_gostshyp_component_api(verror, pressure) result(vcomponent) &
         & bind(C, name=namespace//"new_gostshyp_component")
      !> Error handle
      type(c_ptr), value :: verror
      !> Applied hydrostatic pressure in Hartree/bohr**3
      real(c_double), value :: pressure
      !> New component handle
      type(c_ptr) :: vcomponent
      !> Decoded error handle
      type(vp_error), pointer :: error
      !> Component wrapper
      type(vp_component), pointer :: component
      !> Concrete GOSTSHYP component
      type(solvation_model_component_gostshyp) :: item

      vcomponent = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      allocate (component)
      call new_context(component%ctx, verbosity=0, debug=.false.)
      call new_component_gostshyp(item, real(pressure, wp))
      item%ctx => component%ctx
      allocate (component%ptr, source=item)
      vcomponent = c_loc(component)

   end function new_gostshyp_component_api

   !> Delete a solvation-component handle
   subroutine delete_solvation_component_api(vcomponent) &
         & bind(C, name=namespace//"delete_component")
      !> Component handle
      type(c_ptr), intent(inout), optional :: vcomponent
      !> Decoded component wrapper
      type(vp_component), pointer :: component

      if (.not. present(vcomponent)) return
      if (c_associated(vcomponent)) then
         call c_f_pointer(vcomponent, component)
         if (allocated(component%ptr)) deallocate (component%ptr)
         call component%ctx%delete()
         deallocate (component)
         vcomponent = c_null_ptr
      end if

   end subroutine delete_solvation_component_api

   !> Create a general solvation model around an owned copy of a cavity
   function new_general_solvation_model_api(verror, vcavity, c_debug, c_verbose) result(vmodel)
      !> Error handle
      type(c_ptr), value :: verror
      !> Source cavity handle
      type(c_ptr), value :: vcavity
      !> Debug flag
      logical(c_bool), value :: c_debug
      !> Verbosity level
      integer(c_int), value :: c_verbose
      !> New model handle
      type(c_ptr) :: vmodel
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded cavity wrapper
      type(vp_cavity), pointer :: cavity
      !> Model wrapper
      type(vp_model), pointer :: model
      !> Concrete general model
      type(solvation_model_general) :: general
      !> Constructor error
      type(error_type), allocatable :: model_error

      vmodel = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcavity)) then
         call api_error(error%ptr, "new_model", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcavity, cavity)
      if (.not. associated(cavity%ptr)) then
         call api_error(error%ptr, "new_model", "Cavity is not initialized")
         return
      end if

      allocate (model)
      call new_context(model%ctx, verbosity=int(c_verbose), debug=logical(c_debug))
      call new_model_general(general, cavity%ptr, model%ctx, model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "new_model", model_error%message)
         call model%ctx%delete()
         deallocate (model)
         return
      end if
      allocate (model%ptr, source=general)
      vmodel = c_loc(model)

   end function new_general_solvation_model_api

   !> Append a component to a general model
   subroutine general_model_add_component_api(verror, vmodel, vcomponent) &
         & bind(C, name=namespace//"add_model_component")
      !> Error handle
      type(c_ptr), value :: verror
      !> General-model handle
      type(c_ptr), value :: vmodel
      !> Component handle
      type(c_ptr), value :: vcomponent
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded model wrapper
      type(vp_model), pointer :: model
      !> Decoded component wrapper
      type(vp_component), pointer :: component
      !> Component-addition error
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vmodel) .or. .not. c_associated(vcomponent)) then
         call api_error(error%ptr, "add_model_component", &
                        "Model or component handle is missing")
         return
      end if
      call c_f_pointer(vmodel, model)
      call c_f_pointer(vcomponent, component)
      if (.not. allocated(model%ptr) .or. .not. allocated(component%ptr)) then
         call api_error(error%ptr, "add_model_component", &
                        "Model or component is not initialized")
         return
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         call general%add_component(component%ptr, model_error)
         if (allocated(model_error)) then
            call api_error(error%ptr, "add_model_component", model_error%message)
         end if
      class default
         call api_error(error%ptr, "add_model_component", &
                        "Model is not a general solvation model")
      end select

   end subroutine general_model_add_component_api

   !* ================================================================================= *!
   !*                         Host coupling protocol (C mirror)                         *!
   !* ================================================================================= *!
   !
   ! The C host drives the same exchange as a Fortran host
   !
   ! - one coupling per model (`moist_new_coupling`), a `prepare_*` per phase
   ! - a cursor over the requests with something missing
   !   (`moist_next_coupling_request`), answering each missing output of the
   !   current request by name
   ! - one `get_energy`/`get_response`/`get_gradient` read per phase, handing
   !   the host part back through a response handle
   ! - a cursor over the response items (`moist_next_response_item`), copying
   !   each array of the current item by name
   ! - request kinds, outputs, response items and arrays are addressed by their
   !   scientific names, the same strings the Fortran bindings use
   !
   ! - array reads accept capacities at least as large as the logical grid and
   !   write only its leading entries
   ! - an answer requires the exact `ngrid`; a mismatch is rejected before the
   !   buffer is read and leaves only that output unanswered
   ! - the failure is reported immediately, and the next model read names the
   !   output as missing unless a corrected answer arrives first

   !> Decode a model handle as a general solvation model
   subroutine api_general_model(vmodel, routine, general, error)
      !> Model handle
      type(c_ptr), intent(in) :: vmodel
      !> Calling entry point
      character(len=*), intent(in) :: routine
      !> Pointer to the general model, null on failure
      type(solvation_model_general), pointer, intent(out) :: general
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Decoded model wrapper
      type(vp_model), pointer :: model

      general => null()
      if (.not. c_associated(vmodel)) then
         call api_error(error, routine, "Model handle is missing")
         return
      end if
      call c_f_pointer(vmodel, model)
      if (.not. allocated(model%ptr)) then
         call api_error(error, routine, "Model is not initialized")
         return
      end if
      select type (ptr => model%ptr)
      type is (solvation_model_general)
         general => ptr
      class default
         call api_error(error, routine, "Model is not a general solvation model")
      end select

   end subroutine api_general_model

   !> Decode a coupling handle
   subroutine api_coupling_handle(vcpl, routine, cpl, error)
      !> Coupling handle
      type(c_ptr), intent(in) :: vcpl
      !> Calling entry point
      character(len=*), intent(in) :: routine
      !> Decoded handle, null on failure
      type(vp_coupling), pointer, intent(out) :: cpl
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      cpl => null()
      if (.not. c_associated(vcpl)) then
         call api_error(error, routine, "Coupling handle is missing")
         return
      end if
      call c_f_pointer(vcpl, cpl)

   end subroutine api_coupling_handle

   !> Decode a response handle
   subroutine api_response_handle(vresp, routine, resp, error)
      !> Response handle
      type(c_ptr), intent(in) :: vresp
      !> Calling entry point
      character(len=*), intent(in) :: routine
      !> Decoded handle, null on failure
      type(vp_response), pointer, intent(out) :: resp
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      resp => null()
      if (.not. c_associated(vresp)) then
         call api_error(error, routine, "Response handle is missing")
         return
      end if
      call c_f_pointer(vresp, resp)

   end subroutine api_response_handle

   !> Copy a rank-1 grid array into a caller buffer of the same size
   subroutine api_copy_grid_vector(routine, label, src, c_dst, error)
      !> Calling entry point
      character(len=*), intent(in) :: routine
      !> Name of the array
      character(len=*), intent(in) :: label
      !> Source array, possibly unallocated
      real(wp), allocatable, intent(in) :: src(:)
      !> Destination buffer
      type(c_ptr), intent(in) :: c_dst
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Destination view
      real(c_double), pointer :: dst(:)

      if (.not. allocated(src)) then
         call api_error(error, routine, "'"//label//"' is not available on this coupling")
         return
      end if
      if (.not. c_associated(c_dst)) then
         call api_error(error, routine, "Null array pointer provided for '"//label//"'")
         return
      end if
      call c_f_pointer(c_dst, dst, shape(src))
      dst = real(src, c_double)

   end subroutine api_copy_grid_vector

   !> Copy a rank-2 `(3, ngrid)` grid array into a caller buffer of the same size
   subroutine api_copy_grid_matrix(routine, label, src, c_dst, error)
      !> Calling entry point
      character(len=*), intent(in) :: routine
      !> Name of the array
      character(len=*), intent(in) :: label
      !> Source array, possibly unallocated
      real(wp), allocatable, intent(in) :: src(:, :)
      !> Destination buffer
      type(c_ptr), intent(in) :: c_dst
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Destination view
      real(c_double), pointer :: dst(:, :)

      if (.not. allocated(src)) then
         call api_error(error, routine, "'"//label//"' is not available on this coupling")
         return
      end if
      if (.not. c_associated(c_dst)) then
         call api_error(error, routine, "Null array pointer provided for '"//label//"'")
         return
      end if
      call c_f_pointer(c_dst, dst, shape(src))
      dst = real(src, c_double)

   end subroutine api_copy_grid_matrix

   !> Copy a rank-3 `(3, 3, ngrid)` grid array into an equal-size caller buffer
   subroutine api_copy_grid_tensor(routine, label, src, c_dst, error)
      !> Calling entry point
      character(len=*), intent(in) :: routine
      !> Name of the array
      character(len=*), intent(in) :: label
      !> Source array, possibly unallocated
      real(wp), allocatable, intent(in) :: src(:, :, :)
      !> Destination buffer
      type(c_ptr), intent(in) :: c_dst
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Destination view
      real(c_double), pointer :: dst(:, :, :)

      if (.not. allocated(src)) then
         call api_error(error, routine, "'"//label//"' is not available on this coupling")
         return
      end if
      if (.not. c_associated(c_dst)) then
         call api_error(error, routine, "Null array pointer provided for '"//label//"'")
         return
      end if
      call c_f_pointer(c_dst, dst, shape(src))
      dst = real(src, c_double)

   end subroutine api_copy_grid_tensor

   !> Create the host coupling of a general model
   !>
   !> - runs `new_coupling` on the model: every component declares its requests
   !>   on the updated cavity
   !> - created once per model, survives model updates; the `prepare_*` entry
   !>   points refresh it
   !> - null plus a set error on failure
   function new_coupling_api(verror, vmodel) result(vcpl) &
         & bind(C, name=namespace//"new_coupling")
      !> Error and model handles
      type(c_ptr), value :: verror, vmodel
      !> New coupling handle
      type(c_ptr) :: vcpl
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> General model
      type(solvation_model_general), pointer :: general
      !> New coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Model error
      type(error_type), allocatable :: model_error

      vcpl = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      call api_general_model(vmodel, "new_coupling", general, error%ptr)
      if (allocated(error%ptr)) return

      allocate (cpl)
      call general%new_coupling(cpl%ptr, model_error)
      if (allocated(model_error)) then
         deallocate (cpl)
         call api_error(error%ptr, "new_coupling", model_error%message)
         return
      end if
      cpl%owner => general
      vcpl = c_loc(cpl)

   end function new_coupling_api

   !> Delete a coupling handle and release its model-owned collection
   subroutine delete_coupling_api(vcpl) &
         & bind(C, name=namespace//"delete_coupling")
      !> Coupling handle, null on return
      type(c_ptr), intent(inout), optional :: vcpl
      !> Decoded wrapper
      type(vp_coupling), pointer :: cpl

      if (.not. present(vcpl)) return
      if (c_associated(vcpl)) then
         call c_f_pointer(vcpl, cpl)
         if (associated(cpl%owner)) call cpl%owner%release_coupling(cpl%ptr)
         deallocate (cpl)
         vcpl = c_null_ptr
      end if

   end subroutine delete_coupling_api

   !> Create an empty response handle
   function new_response_api(verror) result(vresp) &
         & bind(C, name=namespace//"new_response")
      !> Error handle
      type(c_ptr), value :: verror
      !> New response handle
      type(c_ptr) :: vresp
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> New response wrapper
      type(vp_response), pointer :: resp

      vresp = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      allocate (resp)
      vresp = c_loc(resp)

   end function new_response_api

   !> Delete a response handle
   subroutine delete_response_api(vresp) &
         & bind(C, name=namespace//"delete_response")
      !> Response handle, null on return
      type(c_ptr), intent(inout), optional :: vresp
      !> Decoded wrapper
      type(vp_response), pointer :: resp

      if (.not. present(vresp)) return
      if (c_associated(vresp)) then
         call c_f_pointer(vresp, resp)
         deallocate (resp)
         vresp = c_null_ptr
      end if

   end subroutine delete_response_api

   !> Decode the handles of a staging entry point
   !>
   !> `cpl` stays null on any failure; `error` stays null for a null error handle
   !>
   !> @param[in]  verror  Error handle
   !> @param[in]  vmodel  Model handle
   !> @param[in]  vcpl    Coupling handle
   !> @param[in]  routine Calling entry point
   !> @param[out] error   Decoded error wrapper
   !> @param[out] general General model
   !> @param[out] cpl     Decoded coupling wrapper
   subroutine api_decode_staging(verror, vmodel, vcpl, routine, error, general, cpl)
      !> Error handle
      type(c_ptr), intent(in) :: verror
      !> Model handle
      type(c_ptr), intent(in) :: vmodel
      !> Coupling handle
      type(c_ptr), intent(in) :: vcpl
      !> Calling entry point
      character(len=*), intent(in) :: routine
      !> Decoded error wrapper
      type(vp_error), pointer, intent(out) :: error
      !> General model, null on failure
      type(solvation_model_general), pointer, intent(out) :: general
      !> Decoded coupling wrapper, null on failure
      type(vp_coupling), pointer, intent(out) :: cpl
      !> Coupling wrapper before validation completes
      type(vp_coupling), pointer :: decoded

      error => null()
      general => null()
      cpl => null()
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call api_general_model(vmodel, routine, general, error%ptr)
      if (allocated(error%ptr)) return
      call api_coupling_handle(vcpl, routine, decoded, error%ptr)
      if (allocated(error%ptr)) return
      cpl => decoded

   end subroutine api_decode_staging

   !> Stage the energy phase: mark every answer stale and restart the host walk
   subroutine general_model_prepare_energy_api(verror, vmodel, vcpl) &
         & bind(C, name=namespace//"prepare_model_energy")
      !> Error, model and coupling handles
      type(c_ptr), value :: verror, vmodel, vcpl
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> General model
      type(solvation_model_general), pointer :: general
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Model error
      type(error_type), allocatable :: model_error

      call api_decode_staging(verror, vmodel, vcpl, "prepare_model_energy", error, general, cpl)
      if (.not. associated(cpl)) return
      call general%prepare_energy(cpl%ptr, model_error)
      if (allocated(model_error)) call api_error(error%ptr, "prepare_model_energy", model_error%message)

   end subroutine general_model_prepare_energy_api

   !> Stage the response phase, retaining valid answers; restarts the host walk
   subroutine general_model_prepare_response_api(verror, vmodel, vcpl) &
         & bind(C, name=namespace//"prepare_model_response")
      !> Error, model and coupling handles
      type(c_ptr), value :: verror, vmodel, vcpl
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> General model
      type(solvation_model_general), pointer :: general
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Model error
      type(error_type), allocatable :: model_error

      call api_decode_staging(verror, vmodel, vcpl, "prepare_model_response", error, general, cpl)
      if (.not. associated(cpl)) return
      call general%prepare_response(cpl%ptr, model_error)
      if (allocated(model_error)) call api_error(error%ptr, "prepare_model_response", model_error%message)

   end subroutine general_model_prepare_response_api

   !> Stage the gradient phase, retaining valid answers; restarts the host walk
   subroutine general_model_prepare_gradient_api(verror, vmodel, vcpl) &
         & bind(C, name=namespace//"prepare_model_gradient")
      !> Error, model and coupling handles
      type(c_ptr), value :: verror, vmodel, vcpl
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> General model
      type(solvation_model_general), pointer :: general
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Model error
      type(error_type), allocatable :: model_error

      call api_decode_staging(verror, vmodel, vcpl, "prepare_model_gradient", error, general, cpl)
      if (.not. associated(cpl)) return
      call general%prepare_gradient(cpl%ptr, model_error)
      if (allocated(model_error)) call api_error(error%ptr, "prepare_model_gradient", model_error%message)

   end subroutine general_model_prepare_gradient_api

   !> Solvation energy of a general model from a staged coupling
   !>
   !> A missing required output of the energy phase, including one whose answer
   !> was rejected, is reported by name and the energy accumulator is unchanged
   subroutine general_model_get_energy_api(verror, vmodel, vcpl, energy) &
         & bind(C, name=namespace//"get_model_energy")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Model handle
      type(c_ptr), value :: vmodel
      !> Coupling handle
      type(c_ptr), value :: vcpl
      !> Energy in Hartree
      real(c_double), intent(inout), optional :: energy
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      !> General model
      type(solvation_model_general), pointer :: general
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Fortran accumulator, added to the caller's only once the call succeeds
      real(wp) :: local
      !> Model error
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(energy)) then
         call api_error(error%ptr, "get_model_energy", "Required pointer 'energy' is missing")
         return
      end if

      call api_general_model(vmodel, "get_model_energy", general, error%ptr)
      if (allocated(error%ptr)) return
      call api_coupling_handle(vcpl, "get_model_energy", cpl, error%ptr)
      if (allocated(error%ptr)) return

      local = 0.0_wp
      call general%get_energy(cpl%ptr, local, model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "get_model_energy", model_error%message)
         return
      end if
      energy = energy + real(local, c_double)

   end subroutine general_model_get_energy_api

   !> Host part of the response phase of a general model from a staged coupling
   !>
   !> - cleared on entry, then filled with the complete host part of the
   !>   response phase: the potential adjoint, the density weights of a field-dependent
   !>   cavity and the GOSTSHYP amplitudes, whichever the model produces
   !> - walk it with `moist_next_response_item` and copy the arrays of each item
   !>   with `moist_get_response_array`
   subroutine general_model_get_response_api(verror, vmodel, vcpl, vresp) &
         & bind(C, name=namespace//"get_model_response")
      !> Error, model, coupling and response handles
      type(c_ptr), value :: verror, vmodel, vcpl, vresp
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> General model
      type(solvation_model_general), pointer :: general
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Decoded response wrapper
      type(vp_response), pointer :: resp
      !> Model error
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      call api_general_model(vmodel, "get_model_response", general, error%ptr)
      if (allocated(error%ptr)) return
      call api_coupling_handle(vcpl, "get_model_response", cpl, error%ptr)
      if (allocated(error%ptr)) return
      call api_response_handle(vresp, "get_model_response", resp, error%ptr)
      if (allocated(error%ptr)) return

      call general%get_response(cpl%ptr, resp%ptr, model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "get_model_response", model_error%message)
      end if

   end subroutine general_model_get_response_api

   !> Nuclear gradient of a general model from a staged coupling
   !>
   !> The response handle is cleared on entry and returns the host part of the
   !> gradient phase (the potential adjoint and GOSTSHYP amplitudes, which the host
   !> contracts with its own geometry derivatives). `gradient` is Fortran
   !> `(3, nat_cap)`; only the leading `nat` columns are written, after the
   !> capacity check
   subroutine general_model_get_gradient_api(verror, vmodel, vcpl, vresp, &
         & nat_cap, c_gradient) &
         & bind(C, name=namespace//"get_model_gradient")
      !> Error, model, coupling and response handles
      type(c_ptr), value :: verror, vmodel, vcpl, vresp
      !> Caller's atom capacity
      integer(c_int), value :: nat_cap
      !> Output gradient, Fortran (3, nat_cap)
      type(c_ptr), value :: c_gradient
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> General model
      type(solvation_model_general), pointer :: general
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Decoded response wrapper
      type(vp_response), pointer :: resp
      !> Output view
      real(c_double), pointer :: gradient(:, :)
      !> Fortran gradient accumulator
      real(wp), allocatable :: local(:, :)
      !> Number of atoms
      integer :: nat
      !> Model error
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      call api_general_model(vmodel, "get_model_gradient", general, error%ptr)
      if (allocated(error%ptr)) return
      call api_coupling_handle(vcpl, "get_model_gradient", cpl, error%ptr)
      if (allocated(error%ptr)) return
      call api_response_handle(vresp, "get_model_gradient", resp, error%ptr)
      if (allocated(error%ptr)) return

      if (.not. general%updated) then
         call api_error(error%ptr, "get_model_gradient", &
                        "General model must be updated first")
         return
      end if
      nat = general%cavity%nsph
      if (.not. c_associated(c_gradient)) then
         call api_error(error%ptr, "get_model_gradient", &
                        "Null gradient pointer provided")
         return
      end if
      if (nat_cap < nat) then
         call api_error(error%ptr, "get_model_gradient", &
                        "Array capacity is too small for 'gradient' - use the atom count of the structure")
         return
      end if

      allocate (local(3, nat), source=0.0_wp)
      call general%get_gradient(cpl%ptr, resp%ptr, local, model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "get_model_gradient", model_error%message)
         return
      end if
      call c_f_pointer(c_gradient, gradient, [3, int(nat_cap)])
      gradient(:, :nat) = gradient(:, :nat) + real(local, c_double)

   end subroutine general_model_get_gradient_api

   !> Install internal model density; update the model before using its results
   subroutine set_model_isodensity_density_api(verror, vmodel, ncart, c_density) &
      bind(C, name=namespace//"set_model_isodensity_density")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Model handle
      type(c_ptr), value, intent(in) :: vmodel
      !> Density dimension
      integer(c_int), value, intent(in) :: ncart
      !> Density buffer
      type(c_ptr), value, intent(in) :: c_density
      !> Error wrapper
      type(vp_error), pointer :: error
      !> Model owner
      type(solvation_model_general), pointer :: general
      !> Native view with reversed C axes
      real(c_double), pointer :: density(:, :)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call api_general_model(vmodel, "set_model_isodensity_density", general, error%ptr)
      if (allocated(error%ptr)) return
      if (ncart < 1 .or. .not. c_associated(c_density)) then
         call api_error(error%ptr, "set_model_isodensity_density", "Invalid density buffer")
         return
      end if
      call c_f_pointer(c_density, density, [int(ncart), int(ncart)])
      call general%set_isodensity_density(real(density, wp), error%ptr)
      call prefix_api_error(error%ptr, "set_model_isodensity_density")
   end subroutine set_model_isodensity_density_api

   !> Advance to the next request with a missing output of the staged phase
   !>
   !> - false ends the pass and rewinds; the next call starts a new pass
   !> - false on any failure too
   function next_coupling_request_api(verror, vcpl) result(more) &
         & bind(C, name=namespace//"next_coupling_request")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Coupling handle
      type(c_ptr), value, intent(in) :: vcpl
      !> Whether a request is now current
      logical(c_bool) :: more
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      more = .false._c_bool
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call api_coupling_handle(vcpl, "next_coupling_request", cpl, error%ptr)
      if (allocated(error%ptr)) return
      more = logical(cpl%ptr%next(), c_bool)
   end function next_coupling_request_api

   !> Scientific name of the current request
   !>
   !> The caller's buffer holds MOIST_NAME_MAX + 1 characters, enough for every
   !> name with its terminator
   subroutine coupling_request_name_api(verror, vcpl, name) &
      bind(C, name=namespace//"get_coupling_request_name")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Coupling handle
      type(c_ptr), value, intent(in) :: vcpl
      !> Name buffer
      character(kind=c_char), intent(inout), optional :: name(*)
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Copy of the current request
      class(coupling_request_type), allocatable :: item
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(name)) then
         call api_error(error%ptr, "get_coupling_request_name", "Output pointer is missing")
         return
      end if
      call api_coupling_handle(vcpl, "get_coupling_request_name", cpl, error%ptr)
      if (allocated(error%ptr)) return
      call current_request(cpl%ptr, item, error%ptr)
      call prefix_api_error(error%ptr, "get_coupling_request_name")
      if (allocated(error%ptr)) return
      call f_c_character(trim(item%name()), name, len_trim(item%name()) + 1)
   end subroutine coupling_request_name_api

   !> Whether an output of the current request is required and still missing
   !>
   !> An output name the request does not declare is not missing
   subroutine coupling_request_missing_api(verror, vcpl, c_name, c_value) &
      bind(C, name=namespace//"get_coupling_request_missing")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Coupling handle
      type(c_ptr), value, intent(in) :: vcpl
      !> NUL-terminated output name
      type(c_ptr), value, intent(in) :: c_name
      !> Required output pointer; unchanged on error
      type(c_ptr), value, intent(in) :: c_value
      !> Caller-owned scalar output
      logical(c_bool), pointer :: output
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Copy of the current request
      class(coupling_request_type), allocatable :: item
      !> Decoded output name
      character(len=:, kind=c_char), allocatable :: name
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. c_associated(c_value)) then
         call api_error(error%ptr, "get_coupling_request_missing", "Output pointer is missing")
         return
      end if
      call api_coupling_handle(vcpl, "get_coupling_request_missing", cpl, error%ptr)
      if (allocated(error%ptr)) return
      call current_request(cpl%ptr, item, error%ptr)
      call prefix_api_error(error%ptr, "get_coupling_request_missing")
      if (allocated(error%ptr)) return
      call api_output_name(c_name, "get_coupling_request_missing", name, error%ptr)
      if (allocated(error%ptr)) return
      call c_f_pointer(c_value, output)
      output = logical(item%is_missing(name), c_bool)
   end subroutine coupling_request_missing_api

   !> Decode a NUL-terminated output or field name argument
   subroutine api_output_name(c_name, routine, name, error, max_len)
      !> NUL-terminated name
      type(c_ptr), intent(in) :: c_name
      !> Calling entry point
      character(len=*), intent(in) :: routine
      !> Decoded name
      character(len=:, kind=c_char), allocatable, intent(out) :: name
      !> Missing or invalid name
      type(error_type), allocatable, intent(out) :: error
      !> Name capacity, output names by default
      integer, intent(in), optional :: max_len
      logical :: truncated
      integer :: cap
      cap = output_name_len
      if (present(max_len)) cap = max_len
      if (.not. c_associated(c_name)) then
         call api_error(error, routine, "Output name is missing")
         return
      end if
      call c_f_character_ptr(c_name, name, cap + 1, truncated)
      if (truncated .or. len(name) == 0) then
         call api_error(error, routine, "Invalid output name")
      end if
   end subroutine api_output_name

   !> Answer one output of the current request
   !>
   !> values is row-major (ngrid, dims...) with the leading extents of the output
   !> and the cavity's grid size; moist reads exactly that many values
   subroutine coupling_answer_api(verror, vcpl, c_name, values) &
      bind(C, name=namespace//"answer_coupling_request")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Coupling handle
      type(c_ptr), value, intent(in) :: vcpl
      !> NUL-terminated output name
      type(c_ptr), value, intent(in) :: c_name
      !> Answer buffer, absent for NULL; assumed size is the one exception to the
      !> c_ptr idiom here: the rows of an output are private to the coupling
      real(c_double), intent(in), optional :: values(*)
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Decoded output name
      character(len=:, kind=c_char), allocatable :: name
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call api_coupling_handle(vcpl, "answer_coupling_request", cpl, error%ptr)
      if (allocated(error%ptr)) return
      call api_output_name(c_name, "answer_coupling_request", name, error%ptr)
      if (allocated(error%ptr)) return
      if (.not. present(values)) then
         call api_error(error%ptr, "answer_coupling_request", "Null array pointer provided for '"//name//"'")
         return
      end if
      call answer_flat(cpl%ptr, name, values, error%ptr)
      call prefix_api_error(error%ptr, "answer_coupling_request")
   end subroutine coupling_answer_api

   !> Copy the Gaussian moment exponents of the current request, ngrid values
   subroutine coupling_get_gaussian_moment_width_api(verror, vcpl, c_values) &
      bind(C, name=namespace//"get_coupling_request_width")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Coupling handle
      type(c_ptr), value, intent(in) :: vcpl
      !> Caller-owned array buffer
      type(c_ptr), value, intent(in) :: c_values
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded coupling wrapper
      type(vp_coupling), pointer :: cpl
      !> Copy of the current request
      class(coupling_request_type), allocatable :: item
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call api_coupling_handle(vcpl, "get_coupling_request_width", cpl, error%ptr)
      if (allocated(error%ptr)) return
      call current_request(cpl%ptr, item, error%ptr)
      call prefix_api_error(error%ptr, "get_coupling_request_width")
      if (allocated(error%ptr)) return
      select type (item)
      type is (gaussian_moment_request_type)
         call api_copy_grid_vector("get_coupling_request_width", "width", item%width, &
                                   & c_values, error%ptr)
      class default
         call api_error(error%ptr, "get_coupling_request_width", &
            & trim(item%name())//" has no input 'width'")
      end select
   end subroutine coupling_get_gaussian_moment_width_api

   !> Advance to the next item of the response
   !>
   !> - false ends the pass and rewinds; the next call starts a new pass
   !> - false on any failure too
   function next_response_item_api(verror, vresp) result(more) &
         & bind(C, name=namespace//"next_response_item")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Response handle
      type(c_ptr), value, intent(in) :: vresp
      !> Whether an item is now current
      logical(c_bool) :: more
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded response wrapper
      type(vp_response), pointer :: resp
      more = .false._c_bool
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call api_response_handle(vresp, "next_response_item", resp, error%ptr)
      if (allocated(error%ptr)) return
      more = logical(resp%ptr%next(), c_bool)
   end function next_response_item_api

   !> Scientific name of the current response item
   !>
   !> The caller's buffer holds MOIST_NAME_MAX + 1 characters, enough for every
   !> name with its terminator
   subroutine response_item_name_api(verror, vresp, name) &
      bind(C, name=namespace//"get_response_item_name")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Response handle
      type(c_ptr), value, intent(in) :: vresp
      !> Name buffer
      character(kind=c_char), intent(inout), optional :: name(*)
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded response wrapper
      type(vp_response), pointer :: resp
      !> Copy of the current item
      class(response_item_type), allocatable :: item
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(name)) then
         call api_error(error%ptr, "get_response_item_name", "Output pointer is missing")
         return
      end if
      call api_response_handle(vresp, "get_response_item_name", resp, error%ptr)
      if (allocated(error%ptr)) return
      call current_response_item(resp%ptr, item, error%ptr)
      call prefix_api_error(error%ptr, "get_response_item_name")
      if (allocated(error%ptr)) return
      call f_c_character(trim(item%name()), name, len_trim(item%name()) + 1)
   end subroutine response_item_name_api

   !> Copy one named array of the current response item
   !>
   !> - C row-major with the grid axis first; writes exactly the array's
   !>   ngrid * dims values into the caller's buffer
   !> - reads the item the cursor stopped at, as an answer writes the current
   !>   request; an array that item does not have is refused by name
   subroutine response_get_api(verror, vresp, c_array, c_values) &
         & bind(C, name=namespace//"get_response_array")
      !> Error handle
      type(c_ptr), value, intent(in) :: verror
      !> Response handle
      type(c_ptr), value, intent(in) :: vresp
      !> NUL-terminated array name
      type(c_ptr), value, intent(in) :: c_array
      !> Output buffer
      type(c_ptr), value, intent(in) :: c_values
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded response wrapper
      type(vp_response), pointer :: resp
      !> Decoded array name
      character(len=:, kind=c_char), allocatable :: array
      !> Copy of the current item
      class(response_item_type), allocatable :: item
      !> Whether the current item has the array
      logical :: known
      !> Public entry point name
      character(len=*), parameter :: routine = "get_response_array"

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      call api_response_handle(vresp, routine, resp, error%ptr)
      if (allocated(error%ptr)) return
      call api_output_name(c_array, routine, array, error%ptr, response_name_len)
      if (allocated(error%ptr)) return
      call current_response_item(resp%ptr, item, error%ptr)
      call prefix_api_error(error%ptr, routine)
      if (allocated(error%ptr)) return
      known = .true.
      select type (item)
      type is (potential_adjoint_response_type)
         select case (array)
         case ("w_phi")
            call api_copy_grid_vector(routine, array, item%w_phi, c_values, error%ptr)
         case default
            known = .false.
         end select
      type is (density_response_type)
         select case (array)
         case ("w_rho")
            call api_copy_grid_vector(routine, array, item%w_rho, c_values, error%ptr)
         case ("w_grad_rho")
            call api_copy_grid_matrix(routine, array, item%w_grad_rho, c_values, error%ptr)
         case ("w_hess_rho")
            call api_copy_grid_tensor(routine, array, item%w_hess_rho, c_values, error%ptr)
         case default
            known = .false.
         end select
      type is (gostshyp_amplitude_response_type)
         select case (array)
         case ("w_overlap")
            call api_copy_grid_vector(routine, array, item%w_overlap, c_values, error%ptr)
         case ("w_normal_deriv")
            call api_copy_grid_vector(routine, array, item%w_normal_deriv, c_values, error%ptr)
         case default
            known = .false.
         end select
      class default
         known = .false.
      end select
      if (.not. known) then
         call api_error(error%ptr, routine, trim(item%name())//" has no array '"//array//"'")
      end if
   end subroutine response_get_api

   !> Report the Cartesian-component layout of an internal isodensity basis
   !>
   !> - query ncart/nshell with NULL arrays, allocate, then call again
   !> - shell offsets are zero-based
   !> - monomial powers define the density-matrix ordering
   subroutine get_isodensity_cart_layout_api(verror, vcav, c_ncart, c_nshell, &
         c_shell_off, c_comp_lx, c_comp_ly, c_comp_lz) &
         & bind(C, name=namespace//"get_isodensity_cart_layout")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Optional output for the number of Cartesian components
      type(c_ptr), value :: c_ncart
      !> Optional output for the number of Gaussian shells
      type(c_ptr), value :: c_nshell
      !> Optional shell-offset buffer, length nshell + 1, zero-based
      type(c_ptr), value :: c_shell_off
      !> Optional x monomial powers, length ncart
      type(c_ptr), value :: c_comp_lx
      !> Optional y monomial powers, length ncart
      type(c_ptr), value :: c_comp_ly
      !> Optional z monomial powers, length ncart
      type(c_ptr), value :: c_comp_lz

      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav
      integer(c_int), pointer :: p_ncart, p_nshell
      integer(c_int), pointer :: shell_off(:), comp_lx(:), comp_ly(:), comp_lz(:)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_isodensity_cart_layout", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)
      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_isodensity_cart_layout", "Cavity is not initialized")
         return
      end if

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         if (.not. allocated(cavity%lsf_model)) then
            call api_error(error%ptr, "get_isodensity_cart_layout", "Cavity has no LSF model")
            return
         end if
         select type (lsf => cavity%lsf_model)
         type is (moist_cavity_drop_lsf_isodensity_internal_type)
            if (c_associated(c_ncart)) then
               call c_f_pointer(c_ncart, p_ncart)
               p_ncart = int(lsf%gto%ncart, c_int)
            end if
            if (c_associated(c_nshell)) then
               call c_f_pointer(c_nshell, p_nshell)
               p_nshell = int(lsf%gto%nshell, c_int)
            end if
            if (c_associated(c_shell_off)) then
               call c_f_pointer(c_shell_off, shell_off, [lsf%gto%nshell + 1])
               shell_off = int(lsf%gto%sh_coff, c_int)
            end if
            if (c_associated(c_comp_lx)) then
               call c_f_pointer(c_comp_lx, comp_lx, [lsf%gto%ncart])
               comp_lx = int(lsf%gto%comp_l(1, :), c_int)
            end if
            if (c_associated(c_comp_ly)) then
               call c_f_pointer(c_comp_ly, comp_ly, [lsf%gto%ncart])
               comp_ly = int(lsf%gto%comp_l(2, :), c_int)
            end if
            if (c_associated(c_comp_lz)) then
               call c_f_pointer(c_comp_lz, comp_lz, [lsf%gto%ncart])
               comp_lz = int(lsf%gto%comp_l(3, :), c_int)
            end if
         class default
            call api_error(error%ptr, "get_isodensity_cart_layout", &
                           "Cavity LSF is not the internal isodensity type")
         end select
      class default
         call api_error(error%ptr, "get_isodensity_cart_layout", &
                        "Cavity is not DROP type")
      end select
   end subroutine get_isodensity_cart_layout_api

   !> Install the cartesian-monomial density matrix for the internal isodensity LSF
   !>
   !> - dcart is the ncart-by-ncart density matrix in moist's cartesian-monomial
   !>   basis (column-major), typically the host's density transformed by the
   !>   cart<-basis map derived from get_isodensity_cart_layout
   !> - call before each cavity build
   subroutine set_isodensity_density_api(verror, vcav, ncart, c_dcart) &
         & bind(C, name=namespace//"set_isodensity_density")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Number of Cartesian components; must match the stored basis
      integer(c_int), value :: ncart
      !> Density matrix in Cartesian-component order, shape (ncart,ncart)
      type(c_ptr), value :: c_dcart

      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav
      real(c_double), pointer :: dcart(:, :)
      type(error_type), allocatable :: cavity_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "set_isodensity_density", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)
      if (.not. cav%owned) then
         call api_error(error%ptr, "set_isodensity_density", "Cannot modify a model-owned cavity")
         return
      end if

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "set_isodensity_density", "Cavity is not initialized")
         return
      end if
      if (.not. c_associated(c_dcart)) then
         call api_error(error%ptr, "set_isodensity_density", "Density pointer is missing")
         return
      end if

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         if (.not. allocated(cavity%lsf_model)) then
            call api_error(error%ptr, "set_isodensity_density", "Cavity has no LSF model")
            return
         end if
         select type (lsf => cavity%lsf_model)
         type is (moist_cavity_drop_lsf_isodensity_internal_type)
            if (int(ncart) /= lsf%gto%ncart) then
               call api_error(error%ptr, "set_isodensity_density", &
                              "Density matrix dimension does not match the basis")
               return
            end if
            call c_f_pointer(c_dcart, dcart, [int(ncart), int(ncart)])
            call lsf%set_density(dcart, cavity_error)
            if (allocated(cavity_error)) then
               call api_error(error%ptr, "set_isodensity_density", cavity_error%message)
               return
            end if
         class default
            call api_error(error%ptr, "set_isodensity_density", &
                           "Cavity LSF is not the internal isodensity type")
         end select
      class default
         call api_error(error%ptr, "set_isodensity_density", "Cavity is not DROP type")
      end select
   end subroutine set_isodensity_density_api

   !> Return the master numerical tolerance configured on a DROP cavity
   subroutine get_drop_cavity_tolerance_api(verror, vcav, c_tolerance) &
         & bind(C, name=namespace//"get_drop_cavity_tolerance")
      !> Opaque error handle supplied by the caller
      type(c_ptr), value :: verror
      !> Decoded API error wrapper
      type(vp_error), pointer :: error
      !> Opaque cavity handle supplied by the caller
      type(c_ptr), value :: vcav
      !> Decoded cavity wrapper
      type(vp_cavity), pointer :: cav
      !> Caller-owned destination for the tolerance
      type(c_ptr), value :: c_tolerance
      !> Decoded destination
      real(c_double), pointer :: tolerance

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_drop_cavity_tolerance", "Cavity handle is missing")
         return
      end if
      if (.not. c_associated(c_tolerance)) then
         call api_error(error%ptr, "get_drop_cavity_tolerance", "Output pointer is missing")
         return
      end if

      call c_f_pointer(vcav, cav)
      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_drop_cavity_tolerance", "Cavity is not initialized")
         return
      end if

      call c_f_pointer(c_tolerance, tolerance)
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         tolerance = real(cavity%param%tolerance, c_double)
      class default
         call api_error(error%ptr, "get_drop_cavity_tolerance", &
                        "Supplied cavity is not a DROP cavity")
      end select
   end subroutine get_drop_cavity_tolerance_api

   !> Assemble A-matrix and compute xi values
   !>
   !> - call before accessing xi values or using the A-matrix
   !> - `ngrid_cap` is the caller's allocated grid capacity, at least the
   !>   cavity's own ngrid
   subroutine assemble_amat_api(verror, vcav, ngrid_cap, amat0, xi) &
         & bind(C, name=namespace//"assemble_amat")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Allocated grid-point capacity
      integer(c_int), value :: ngrid_cap
      !> PCM interaction matrix
      real(c_double), intent(inout), optional :: amat0(ngrid_cap, ngrid_cap)
      !> Gaussian inverse widths
      real(c_double), intent(inout), optional :: xi(ngrid_cap)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav
      real(wp), allocatable :: amat0_local(:, :)
      type(error_type), allocatable :: cavity_error
      integer :: ngrid

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(amat0)) then
         call api_error(error%ptr, "assemble_amat", "Required pointer 'amat0' is missing")
         return
      end if
      if (.not. present(xi)) then
         call api_error(error%ptr, "assemble_amat", "Required pointer 'xi' is missing")
         return
      end if
      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "assemble_amat", "DROP cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "assemble_amat", "DROP cavity is not initialized")
         return
      end if

      associate (cavity => cav%ptr)
         ngrid = cavity%ngrid
         if (ngrid_cap < ngrid) then
            call api_error(error%ptr, "assemble_amat", &
                           "Array capacity is too small - use ngrid value from get_cavity_sizes")
            return
         end if
         if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f) &
             .or. .not. allocated(cavity%xyz)) then
            call api_error(error%ptr, "assemble_amat", &
                           "Cavity does not expose Gaussian surface widths")
            return
         end if

         ! Generic Gaussian surface-charge interaction matrix
         allocate (amat0_local(ngrid, ngrid))
         call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, amat0_local, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "assemble_amat", cavity_error%message)
            return
         end if

         ! Copy results to output arrays
         amat0(:ngrid, :ngrid) = amat0_local(:, :)
         xi(:cavity%ngrid) = cavity%xi0
      end associate

   end subroutine assemble_amat_api

   !> Return Gaussian widths and switching factors for a Gaussian PCM cavity
   subroutine get_cavity_gaussian_api(verror, vcav, ngrid_cap, c_xi, c_f) &
         & bind(C, name=namespace//"get_cavity_gaussian")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Cavity handle
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      !> Allocated grid capacity; must cover the logical grid
      integer(c_int), value :: ngrid_cap
      !> Output Gaussian widths, length ngrid_cap
      type(c_ptr), value :: c_xi
      !> Output switching factors, length ngrid_cap
      type(c_ptr), value :: c_f
      real(c_double), pointer :: xi(:), f(:)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_cavity_gaussian", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_cavity_gaussian", "Cavity is not initialized")
         return
      end if
      if (ngrid_cap < cav%ptr%ngrid) then
         call api_error(error%ptr, "get_cavity_gaussian", &
                        "Array capacity is too small - use ngrid value from get_cavity_sizes")
         return
      end if
      if (.not. c_associated(c_xi) .or. .not. c_associated(c_f)) then
         call api_error(error%ptr, "get_cavity_gaussian", "Null array pointer provided")
         return
      end if
      if (.not. allocated(cav%ptr%xi0) .or. .not. allocated(cav%ptr%f)) then
         call api_error(error%ptr, "get_cavity_gaussian", &
                        "Cavity does not provide a Gaussian PCM surface")
         return
      end if
      if (size(cav%ptr%xi0) /= cav%ptr%ngrid .or. size(cav%ptr%f) /= cav%ptr%ngrid) then
         call api_error(error%ptr, "get_cavity_gaussian", &
                        "Gaussian PCM surface has inconsistent dimensions")
         return
      end if

      call c_f_pointer(c_xi, xi, [int(ngrid_cap)])
      call c_f_pointer(c_f, f, [int(ngrid_cap)])
      xi(:cav%ptr%ngrid) = cav%ptr%xi0
      f(:cav%ptr%ngrid) = cav%ptr%f

   end subroutine get_cavity_gaussian_api

   !* ================================================================================= *!
   !*                            Generic cavity API (Tier 1)                            *!
   !* ================================================================================= *!

   !> Update a cavity of any concrete type
   subroutine update_cavity_api(verror, vcav, vmol) &
         & bind(C, name=namespace//"update_cavity")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Cavity handle
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      !> Molecular structure handle
      type(c_ptr), value :: vmol
      type(vp_structure), pointer :: mol
      type(error_type), allocatable :: cavity_error
      integer :: nat

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "update_cavity", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "update_cavity", "Cavity is not initialized")
         return
      end if
      if (.not. cav%owned) then
         call api_error(error%ptr, "update_cavity", &
                        "Cannot update a borrowed cavity handle")
         return
      end if

      if (.not. c_associated(vmol)) then
         call api_error(error%ptr, "update_cavity", "Molecular structure data is missing")
         return
      end if
      call c_f_pointer(vmol, mol)

      nat = mol%ptr%nat
      if (nat <= 0) then
         call api_error(error%ptr, "update_cavity", "Invalid number of atoms")
         return
      end if

      ! Dispatch through the deferred update procedure
      call cav%ptr%update(mol%ptr, error=cavity_error)
      if (allocated(cavity_error)) then
         call api_error(error%ptr, "update_cavity", cavity_error%message)
         return
      end if

   end subroutine update_cavity_api

   !> Report grid and sphere counts of a cavity of any concrete type
   subroutine get_cavity_sizes_api(verror, vcav, ngrid, nsph) &
         & bind(C, name=namespace//"get_cavity_sizes")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Number of surface grid points
      integer(c_int), intent(inout), optional :: ngrid
      !> Number of atomic spheres
      integer(c_int), intent(inout), optional :: nsph
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      !> Result staged until the operation succeeds
      integer(c_int) :: local_ngrid
      !> Result staged until the operation succeeds
      integer(c_int) :: local_nsph
      type(vp_cavity), pointer :: cav

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(ngrid)) then
         call api_error(error%ptr, "get_cavity_sizes", "Required pointer 'ngrid' is missing")
         return
      end if
      if (.not. present(nsph)) then
         call api_error(error%ptr, "get_cavity_sizes", "Required pointer 'nsph' is missing")
         return
      end if
      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_cavity_sizes", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_cavity_sizes", "Cavity is not initialized")
         return
      end if

      if (.not. allocated(cav%ptr%radii)) then
         call api_error(error%ptr, "get_cavity_sizes", "Cavity is not built yet")
         return
      end if

      ! Base cavity_type fields
      local_ngrid = cav%ptr%ngrid
      local_nsph = size(cav%ptr%radii)
      if (allocated(error%ptr)) return
      ngrid = local_ngrid
      nsph = local_nsph

   end subroutine get_cavity_sizes_api

   !> Report cavity results common to any concrete cavity type
   !>
   !> - only fields from the base cavity_type
   !> - the caller passes the capacities it allocated the arrays with, at least
   !>   the cavity's own ngrid/nsph, otherwise nothing is written
   subroutine get_cavity_results_api(verror, vcav, ngrid_cap, nsph_cap, &
         & area, volume, ngrid, nsph, &
         & xyz, a, owner, converged, radii, asph) &
         & bind(C, name=namespace//"get_cavity_results")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Allocated grid-point capacity
      integer(c_int), value :: ngrid_cap
      !> Allocated sphere capacity
      integer(c_int), value :: nsph_cap
      !> Total cavity area in Bohr squared
      real(c_double), intent(inout), optional :: area
      !> Enclosed cavity volume in Bohr cubed
      real(c_double), intent(inout), optional :: volume
      !> Number of surface grid points
      integer(c_int), intent(inout), optional :: ngrid
      !> Number of atomic spheres
      integer(c_int), intent(inout), optional :: nsph
      !> Surface coordinates, shape (3,ngrid_cap)
      real(c_double), intent(inout), optional :: xyz(3, ngrid_cap)
      !> PCM matrix in column-major order
      real(c_double), intent(inout), optional :: a(ngrid_cap)
      !> Owning atom index for each grid point
      integer(c_int), intent(inout), optional :: owner(ngrid_cap)
      !> Whether the cavity projection converged
      logical(c_bool), intent(inout), optional :: converged(ngrid_cap)
      !> Atomic radii in Bohr
      real(c_double), intent(inout), optional :: radii(nsph_cap)
      !> Per-atom surface areas
      real(c_double), intent(inout), optional :: asph(nsph_cap)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      !> Result staged until the operation succeeds
      real(c_double) :: local_area
      !> Result staged until the operation succeeds
      real(c_double) :: local_volume
      !> Result staged until the operation succeeds
      integer(c_int) :: local_ngrid
      !> Result staged until the operation succeeds
      integer(c_int) :: local_nsph
      type(vp_cavity), pointer :: cav

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(area)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'area' is missing")
         return
      end if
      if (.not. present(volume)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'volume' is missing")
         return
      end if
      if (.not. present(ngrid)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'ngrid' is missing")
         return
      end if
      if (.not. present(nsph)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'nsph' is missing")
         return
      end if
      if (.not. present(xyz)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'xyz' is missing")
         return
      end if
      if (.not. present(a)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'a' is missing")
         return
      end if
      if (.not. present(owner)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'owner' is missing")
         return
      end if
      if (.not. present(converged)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'converged' is missing")
         return
      end if
      if (.not. present(radii)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'radii' is missing")
         return
      end if
      if (.not. present(asph)) then
         call api_error(error%ptr, "get_cavity_results", "Required pointer 'asph' is missing")
         return
      end if
      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_cavity_results", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_cavity_results", "Cavity is not initialized")
         return
      end if

      if (.not. allocated(cav%ptr%total_area)) then
         call api_error(error%ptr, "get_cavity_results", "Cavity is not built yet")
         return
      end if

      if (ngrid_cap < cav%ptr%ngrid .or. nsph_cap < size(cav%ptr%radii)) then
         call api_error(error%ptr, "get_cavity_results", &
                        "Array capacity is too small - use get_cavity_sizes")
         return
      end if

      ! Get scalar values from base cavity_type
      local_area = cav%ptr%total_area
      local_volume = cav%ptr%total_volume
      local_ngrid = cav%ptr%ngrid
      local_nsph = size(cav%ptr%radii)

      ! Get grid point arrays from base cavity_type
      xyz(:, :local_ngrid) = cav%ptr%xyz(:, :)
      a(:local_ngrid) = cav%ptr%a
      ! Convert from Fortran 1-based to C 0-based indexing
      owner(:local_ngrid) = cav%ptr%owner - 1
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         converged(:local_ngrid) = cavity%converged
      class default
         converged(:local_ngrid) = .true.
      end select

      ! Get per-sphere arrays from base cavity_type
      radii(:local_nsph) = cav%ptr%radii
      asph(:local_nsph) = cav%ptr%asph
      if (allocated(error%ptr)) return
      area = local_area
      volume = local_volume
      ngrid = local_ngrid
      nsph = local_nsph

   end subroutine get_cavity_results_api

   !* ================================================================================= *!
   !*                            Named cavity fields (Tier 2)                           *!
   !* ================================================================================= *!

   !> Number of named result fields the cavity currently holds
   !>
   !> - built from the cavity's own declarations, so it grows with the cavity
   !>   type and shrinks when optional properties were not requested
   !> - fields keep the names the cavity uses internally
   subroutine get_cavity_field_count_api(verror, vcav, nfield) &
         & bind(C, name=namespace//"get_cavity_field_count")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Number of available fields
      integer(c_int), intent(inout), optional :: nfield
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      !> Result staged until the operation succeeds
      integer(c_int) :: local_nfield
      !> Fortran cavity pointer
      type(vp_cavity), pointer :: cav
      type(cavity_field_query_type) :: query

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(nfield)) then
         call api_error(error%ptr, "get_cavity_field_count", "Required pointer 'nfield' is missing")
         return
      end if
      local_nfield = 0
      if (.not. resolve_field_cavity(verror, vcav, "get_cavity_field_count", error, cav)) return

      call query%enumerate()
      call cav%ptr%list_fields(query)
      local_nfield = query%nfield
      if (allocated(error%ptr)) return
      nfield = local_nfield

   end subroutine get_cavity_field_count_api

   !> Describe one readable field by position
   !>
   !> - `dtype` one of the MOIST_FIELD_* tags
   !> - `rank` 0 for a scalar
   !> - `dims` extents in C row-major order, slowest-varying first
   !> - `count` number of elements a fetch writes
   !> - name buffer holds MOIST_FIELD_NAME_MAX + 1 characters
   subroutine get_cavity_field_info_api(verror, vcav, ifield, name, &
         & dtype, rank, dims, count) &
         & bind(C, name=namespace//"get_cavity_field_info")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Zero-based field index
      integer(c_int), value :: ifield
      !> Field name buffer
      character(kind=c_char), intent(inout), optional :: name(*)
      !> Field scalar type selector
      integer(c_int), intent(inout), optional :: dtype
      !> Number of field dimensions
      integer(c_int), intent(inout), optional :: rank
      !> Field dimensions
      integer(c_int), intent(inout), optional :: dims(cavity_field_max_rank)
      !> Number of available entries
      integer(c_int), intent(inout), optional :: count
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      !> Result staged until the operation succeeds
      integer(c_int) :: local_dtype
      !> Result staged until the operation succeeds
      integer(c_int) :: local_rank
      !> Result staged until the operation succeeds
      integer(c_int) :: local_count
      !> Fortran cavity pointer
      type(vp_cavity), pointer :: cav
      type(cavity_field_query_type) :: query

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(name)) then
         call api_error(error%ptr, "get_cavity_field_info", "Required pointer 'name' is missing")
         return
      end if
      if (.not. present(dtype)) then
         call api_error(error%ptr, "get_cavity_field_info", "Required pointer 'dtype' is missing")
         return
      end if
      if (.not. present(rank)) then
         call api_error(error%ptr, "get_cavity_field_info", "Required pointer 'rank' is missing")
         return
      end if
      if (.not. present(dims)) then
         call api_error(error%ptr, "get_cavity_field_info", "Required pointer 'dims' is missing")
         return
      end if
      if (.not. present(count)) then
         call api_error(error%ptr, "get_cavity_field_info", "Required pointer 'count' is missing")
         return
      end if
      local_dtype = 0
      local_rank = 0
      local_count = 0
      if (.not. resolve_field_cavity(verror, vcav, "get_cavity_field_info", error, cav)) return

      call query%enumerate()
      call cav%ptr%list_fields(query)

      if (ifield < 0 .or. ifield >= query%nfield) then
         call api_error(error%ptr, "get_cavity_field_info", &
            & "Field index out of range - use the count from get_cavity_field_count")
         return
      end if

      associate (info => query%info(ifield + 1))
         if (len(info%name) > max_field_name_len) then
            call api_error(error%ptr, "get_cavity_field_info", "Field name exceeds MOIST_FIELD_NAME_MAX")
            return
         end if
         call f_c_character(info%name, name, len(info%name) + 1)
         local_dtype = info%dtype
         local_rank = info%rank
         dims = 1_c_int
         dims(:local_rank) = info%dims(local_rank:1:-1)
         local_count = info%count()
      end associate
      if (allocated(error%ptr)) return
      dtype = local_dtype
      rank = local_rank
      count = local_count

   end subroutine get_cavity_field_info_api

   !> Copy a field description or query its full length
   subroutine get_cavity_field_about_api(verror, vcav, cname, about, capacity, length) &
         & bind(C, name="moist_get_cavity_field_about")
      !> Required diagnostic handle
      type(c_ptr), value, intent(in) :: verror
      !> Cavity handle
      type(c_ptr), value, intent(in) :: vcav
      !> NUL-terminated field name
      type(c_ptr), value, intent(in) :: cname
      !> Caller-owned description buffer
      character(kind=c_char), intent(inout), optional :: about(*)
      !> Description buffer capacity including terminator
      integer(c_size_t), value, intent(in) :: capacity
      !> Full text length excluding the terminator
      integer(c_size_t), intent(inout), optional :: length
      !> Decoded diagnostic handle
      type(vp_error), pointer :: error
      !> Decoded cavity handle
      type(vp_cavity), pointer :: cav
      !> Field metadata and payload
      type(cavity_field_query_type) :: query

      if (.not. valid_string_output(verror, about, capacity, length, "get_cavity_field_about", error)) return
      if (.not. resolve_field_cavity(verror, vcav, "get_cavity_field_about", error, cav)) return
      if (.not. fetch_cavity_field(error, cav, cname, "get_cavity_field_about", query)) return
      call copy_string_output(query%hit%about, about, capacity, length)
   end subroutine get_cavity_field_about_api

   !> Read a real-valued field by name
   !>
   !> - buffer filled in Fortran order for rank-2 fields
   !> - receives the `count` elements `get_cavity_field_info` reports
   subroutine get_cavity_field_real_api(verror, vcav, cname, values) &
         & bind(C, name=namespace//"get_cavity_field_real")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> NUL-terminated field name
      type(c_ptr), value :: cname
      !> Packed field values
      real(c_double), intent(inout), optional :: values(*)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      !> Fortran cavity pointer
      type(vp_cavity), pointer :: cav
      type(cavity_field_query_type) :: query

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(values)) then
         call api_error(error%ptr, "get_cavity_field_real", "Required pointer 'values' is missing")
         return
      end if
      if (.not. resolve_field_cavity(verror, vcav, "get_cavity_field_real", error, cav)) return
      if (.not. fetch_cavity_field(error, cav, cname, "get_cavity_field_real", query)) return
      if (.not. check_field_payload(error, query, cavity_field_real, &
         & "get_cavity_field_real")) return

      values(:size(query%rvals)) = query%rvals

   end subroutine get_cavity_field_real_api

   !> Read an integer-valued field by name
   !>
   !> Fields that report a sphere index are handed out 0-based, matching
   !> get_cavity_results; ids such as `numbering` and `anchor_id` are passed
   !> through as the cavity stores them
   subroutine get_cavity_field_int_api(verror, vcav, cname, values) &
         & bind(C, name=namespace//"get_cavity_field_int")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> NUL-terminated field name
      type(c_ptr), value :: cname
      !> Packed field values
      integer(c_int), intent(inout), optional :: values(*)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      !> Fortran cavity pointer
      type(vp_cavity), pointer :: cav
      type(cavity_field_query_type) :: query

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(values)) then
         call api_error(error%ptr, "get_cavity_field_int", "Required pointer 'values' is missing")
         return
      end if
      if (.not. resolve_field_cavity(verror, vcav, "get_cavity_field_int", error, cav)) return
      if (.not. fetch_cavity_field(error, cav, cname, "get_cavity_field_int", query)) return
      if (.not. check_field_payload(error, query, cavity_field_int, &
         & "get_cavity_field_int")) return

      values(:size(query%ivals)) = query%ivals

   end subroutine get_cavity_field_int_api

   !> Read a logical-valued field by name
   subroutine get_cavity_field_bool_api(verror, vcav, cname, values) &
         & bind(C, name=namespace//"get_cavity_field_bool")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> NUL-terminated field name
      type(c_ptr), value :: cname
      !> Packed field values
      logical(c_bool), intent(inout), optional :: values(*)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      !> Fortran cavity pointer
      type(vp_cavity), pointer :: cav
      type(cavity_field_query_type) :: query

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(values)) then
         call api_error(error%ptr, "get_cavity_field_bool", "Required pointer 'values' is missing")
         return
      end if
      if (.not. resolve_field_cavity(verror, vcav, "get_cavity_field_bool", error, cav)) return
      if (.not. fetch_cavity_field(error, cav, cname, "get_cavity_field_bool", query)) return
      if (.not. check_field_payload(error, query, cavity_field_bool, &
         & "get_cavity_field_bool")) return

      values(:size(query%lvals)) = logical(query%lvals, c_bool)

   end subroutine get_cavity_field_bool_api

   !> Resolve the error and cavity handles shared by the field entry points
   !>
   !> @return             Whether both handles resolved to a usable cavity
   logical function resolve_field_cavity(verror, vcav, origin, error, cav) result(ok)
      !> Error handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Entry point name used in error messages
      character(len=*), intent(in) :: origin
      !> Fortran error pointer
      type(vp_error), pointer, intent(out) :: error
      !> Fortran cavity pointer
      type(vp_cavity), pointer, intent(out) :: cav

      ok = .false.
      nullify (error)
      nullify (cav)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, origin, "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, origin, "Cavity is not initialized")
         return
      end if

      ok = .true.

   end function resolve_field_cavity

   !> Look one named field up on a cavity
   !>
   !> - an undeclared name is an error, and so is a field whose optional
   !>   property was never requested
   !> - the cavity does not declare an array it has not computed
   !>
   !> @return              Whether the field was found
   logical function fetch_cavity_field(error, cav, cname, origin, query) result(ok)
      !> Fortran error pointer
      type(vp_error), pointer, intent(in) :: error
      !> Fortran cavity pointer
      type(vp_cavity), pointer, intent(in) :: cav
      !> Field name as a C string
      type(c_ptr), value :: cname
      !> Entry point name used in error messages
      character(len=*), intent(in) :: origin
      !> Query walker
      type(cavity_field_query_type), intent(inout) :: query

      character(len=:, kind=c_char), allocatable :: name

      ok = .false.

      if (.not. c_associated(cname)) then
         call api_error(error%ptr, origin, "Field name is missing")
         return
      end if
      call c_f_character_ptr(cname, name, max_field_name_len)

      if (len(name) == 0) then
         call api_error(error%ptr, origin, "Field name is empty")
         return
      end if

      call query%fetch(name)
      call cav%ptr%list_fields(query)

      if (.not. query%found) then
         call api_error(error%ptr, origin, &
            & "Cavity has no field '"//name//"' - it is either unknown or was not computed; "// &
            & "enumerate the available fields with get_cavity_field_count/get_cavity_field_info")
         return
      end if

      ok = .true.

   end function fetch_cavity_field

   !> Check that a fetched field matches the requested element type and fits
   !>
   !> @return               Whether the payload may be copied out
   logical function check_field_payload(error, query, dtype, origin) result(ok)
      !> Fortran error pointer
      type(vp_error), pointer, intent(in) :: error
      !> Query walker holding the fetched payload
      type(cavity_field_query_type), intent(in) :: query
      !> Element type the caller asked for
      integer, intent(in) :: dtype
      !> Entry point name used in error messages
      character(len=*), intent(in) :: origin

      ok = .false.

      if (query%hit%dtype /= dtype) then
         call api_error(error%ptr, origin, &
            & "Field '"//query%hit%name//"' has a different element type - "// &
            & "read the type tag from get_cavity_field_info")
         return
      end if

      ok = .true.

   end function check_field_payload

   !* ================================================================================= *!
   !*                       Cavity and A-matrix gradients (Tier 3)                      *!
   !* ================================================================================= *!

   !> Compute cavity gradient w.r.t. nuclear coordinates
   !>
   !> - call after update_cavity and before get_cavity_gradient
   subroutine compute_cavity_gradient_api(verror, vcav) &
         & bind(C, name=namespace//"compute_cavity_gradient")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Cavity handle
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(error_type), allocatable :: cavity_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "compute_cavity_gradient", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "compute_cavity_gradient", "Cavity is not initialized")
         return
      end if

      if (.not. allocated(cav%ptr%total_area)) then
         call api_error(error%ptr, "compute_cavity_gradient", "Cavity is not built yet - call update_cavity first")
         return
      end if

      ! Enable optional gradient arrays required by get_cavity_gradient
      select type (c => cav%ptr)
      type is (cavity_type_drop)
         c%request%r_iI = .true.
         c%request%rho = .true.
      end select

      ! Call the deferred get_gradient procedure
      call cav%ptr%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call api_error(error%ptr, "compute_cavity_gradient", cavity_error%message)
         return
      end if

   end subroutine compute_cavity_gradient_api

   !> Compute anchor-only nuclear derivatives (callback/isodensity LSF)
   !>
   !> - call after update_cavity and before the *_rA contractions
   !> - restricts each grid point's nuclear coupling to its owner atom's rigid
   !>   anchor motion (the field's nuclear derivatives vanish for callback LSFs)
   subroutine compute_anchor_gradient_api(verror, vcav) &
         & bind(C, name=namespace//"compute_anchor_gradient")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Cavity handle
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(error_type), allocatable :: cavity_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "compute_anchor_gradient", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "compute_anchor_gradient", "Cavity is not initialized")
         return
      end if

      if (.not. allocated(cav%ptr%total_area)) then
         call api_error(error%ptr, "compute_anchor_gradient", "Cavity is not built yet - call update_cavity first")
         return
      end if

      select type (c => cav%ptr)
      type is (cavity_type_drop)
         call c%compute_anchor_gradient(cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "compute_anchor_gradient", cavity_error%message)
            return
         end if
      class default
         call api_error(error%ptr, "compute_anchor_gradient", &
                        "Cavity is not DROP type - anchor gradient only supports DROP cavities")
      end select

   end subroutine compute_anchor_gradient_api

   !> Get the anchor-channel nuclear derivatives produced by the anchor pass
   !>
   !> - call compute_anchor_gradient (or compute_cavity_gradient) first
   !> - arrays in Fortran shape notation (column-major: the leftmost index is
   !>   contiguous); a C caller passes flat buffers of the same total size and
   !>   indexes element (i1,...,in) at i1 + d1*(i2 + d2*(...)), all 0-based:
   !>   xyz1_rA(3, 3, nsph, ngrid)  - d(r_i)_j / d(R_A)_alpha  (j, alpha, A, grid)
   !>   xi1_rA(3, nsph, ngrid)      - d(xi_i)  / d(R_A)_alpha   (alpha, A, grid)
   !>   a_i1_rA(3, nsph, ngrid)     - d(a_i)   / d(R_A)_alpha   (alpha, A, grid)
   !>   v_i1_rA(3, nsph, ngrid)     - d(v_i)   / d(R_A)_alpha   (alpha, A, grid)
   !>   A_tot1_rA(3, nsph)          - d(total area)   / d(R_A)_alpha  (alpha, A)
   !>   V_tot1_rA(3, nsph)          - d(total volume) / d(R_A)_alpha  (alpha, A)
   !> - the per-point area/volume elements (a_i1_rA, v_i1_rA) are the un-summed
   !>   counterparts of A_tot1_rA/V_tot1_rA; summing over the grid recovers the
   !>   totals
   !> - the area carries a switching-function dependence (a_i ~ f_i / xi_i^2),
   !>   so a_i1_rA is NOT recoverable from xi1_rA alone
   !> - used by geometric surface functionals such as GOSTSHYP for their area route
   subroutine get_anchor_gradient_api(verror, vcav, nsph_cap, ngrid_cap, xyz1_rA, xi1_rA, &
         & a_i1_rA, v_i1_rA, A_tot1_rA, V_tot1_rA) &
         & bind(C, name=namespace//"get_anchor_gradient")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Allocated sphere capacity
      integer(c_int), value :: nsph_cap
      !> Allocated grid-point capacity
      integer(c_int), value :: ngrid_cap
      !> Surface-position derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: xyz1_rA(3, 3, nsph_cap, ngrid_cap)
      !> Gaussian-width derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: xi1_rA(3, nsph_cap, ngrid_cap)
      !> Per-point area derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: a_i1_rA(3, nsph_cap, ngrid_cap)
      !> Per-point volume derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: v_i1_rA(3, nsph_cap, ngrid_cap)
      !> Total area derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: A_tot1_rA(3, nsph_cap)
      !> Total volume derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: V_tot1_rA(3, nsph_cap)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(xyz1_rA)) then
         call api_error(error%ptr, "get_anchor_gradient", "Required pointer 'xyz1_rA' is missing")
         return
      end if
      if (.not. present(xi1_rA)) then
         call api_error(error%ptr, "get_anchor_gradient", "Required pointer 'xi1_rA' is missing")
         return
      end if
      if (.not. present(a_i1_rA)) then
         call api_error(error%ptr, "get_anchor_gradient", "Required pointer 'a_i1_rA' is missing")
         return
      end if
      if (.not. present(v_i1_rA)) then
         call api_error(error%ptr, "get_anchor_gradient", "Required pointer 'v_i1_rA' is missing")
         return
      end if
      if (.not. present(A_tot1_rA)) then
         call api_error(error%ptr, "get_anchor_gradient", "Required pointer 'A_tot1_rA' is missing")
         return
      end if
      if (.not. present(V_tot1_rA)) then
         call api_error(error%ptr, "get_anchor_gradient", "Required pointer 'V_tot1_rA' is missing")
         return
      end if
      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_anchor_gradient", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_anchor_gradient", "Cavity is not initialized")
         return
      end if

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         if (.not. allocated(cavity%xyz1_rA) .or. .not. allocated(cavity%xi1_rA) &
             .or. .not. allocated(cavity%a_i1_rA) .or. .not. allocated(cavity%v1_rA) &
             .or. .not. allocated(cavity%A_tot1_rA) .or. .not. allocated(cavity%V_tot1_rA)) then
            call api_error(error%ptr, "get_anchor_gradient", &
                           "Gradient not computed - call compute_anchor_gradient first")
            return
         end if
         if (nsph_cap < cavity%nsph .or. ngrid_cap < cavity%ngrid) then
            call api_error(error%ptr, "get_anchor_gradient", &
                           "Array capacity is too small - use get_cavity_sizes")
            return
         end if
         associate (nsph => cavity%nsph, ngrid => cavity%ngrid)
            xyz1_rA(:, :, :nsph, :ngrid) = cavity%xyz1_rA(:, :, :, :)
            xi1_rA(:, :nsph, :ngrid) = cavity%xi1_rA(:, :, :)
            a_i1_rA(:, :nsph, :ngrid) = cavity%a_i1_rA(:, :, :)
            v_i1_rA(:, :nsph, :ngrid) = cavity%v1_rA(:, :, :)
            A_tot1_rA(:, :nsph) = cavity%A_tot1_rA(:, :)
            V_tot1_rA(:, :nsph) = cavity%V_tot1_rA(:, :)
         end associate
      class default
         call api_error(error%ptr, "get_anchor_gradient", &
                        "Cavity is not DROP type - anchor gradient only supports DROP cavities")
      end select

   end subroutine get_anchor_gradient_api

   !> Get cavity gradient arrays (DROP-specific)
   !>
   !> - call compute_cavity_gradient first
   !> - arrays in Fortran shape notation (column-major: the leftmost index is
   !>   contiguous); a C caller passes flat buffers of the same total size:
   !>   A_tot1_rA(3, nsph)           - gradient of total area
   !>   V_tot1_rA(3, nsph)           - gradient of total volume
   !>   asph1_rA(3, nsph, nsph)      - gradient of per-sphere areas
   !>   vsph1_rA(3, nsph, nsph)      - gradient of per-sphere volumes
   !>   xyz1_rA(3, 3, nsph, ngrid) - grid point position derivatives (j, alpha, A, grid)
   !>   r_iI1_rA(3, nsph, ngrid)     - gradient of grid-owner distances
   !>   rho1_rA(3, nsph, ngrid)      - gradient of rho values
   subroutine get_cavity_gradient_api(verror, vcav, nsph_cap, ngrid_cap, &
         & A_tot1_rA, V_tot1_rA, asph1_rA, vsph1_rA, &
         & xyz1_rA, r_iI1_rA, rho1_rA) &
         & bind(C, name=namespace//"get_cavity_gradient")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Allocated sphere capacity
      integer(c_int), value :: nsph_cap
      !> Allocated grid-point capacity
      integer(c_int), value :: ngrid_cap
      !> Total area derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: A_tot1_rA(3, nsph_cap)
      !> Total volume derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: V_tot1_rA(3, nsph_cap)
      !> Atomic area derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: asph1_rA(3, nsph_cap, nsph_cap)
      !> Atomic volume derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: vsph1_rA(3, nsph_cap, nsph_cap)
      !> Surface-position derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: xyz1_rA(3, 3, nsph_cap, ngrid_cap)
      !> Point-to-atom displacement derivatives
      real(c_double), intent(inout), optional :: r_iI1_rA(3, nsph_cap, ngrid_cap)
      !> Radial-distance derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: rho1_rA(3, nsph_cap, ngrid_cap)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav
      integer :: nsph, ngrid

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(A_tot1_rA)) then
         call api_error(error%ptr, "get_cavity_gradient", "Required pointer 'A_tot1_rA' is missing")
         return
      end if
      if (.not. present(V_tot1_rA)) then
         call api_error(error%ptr, "get_cavity_gradient", "Required pointer 'V_tot1_rA' is missing")
         return
      end if
      if (.not. present(asph1_rA)) then
         call api_error(error%ptr, "get_cavity_gradient", "Required pointer 'asph1_rA' is missing")
         return
      end if
      if (.not. present(vsph1_rA)) then
         call api_error(error%ptr, "get_cavity_gradient", "Required pointer 'vsph1_rA' is missing")
         return
      end if
      if (.not. present(xyz1_rA)) then
         call api_error(error%ptr, "get_cavity_gradient", "Required pointer 'xyz1_rA' is missing")
         return
      end if
      if (.not. present(r_iI1_rA)) then
         call api_error(error%ptr, "get_cavity_gradient", "Required pointer 'r_iI1_rA' is missing")
         return
      end if
      if (.not. present(rho1_rA)) then
         call api_error(error%ptr, "get_cavity_gradient", "Required pointer 'rho1_rA' is missing")
         return
      end if
      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_cavity_gradient", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_cavity_gradient", "Cavity is not initialized")
         return
      end if

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         ! Check if gradient was computed
         if (.not. allocated(cavity%A_tot1_rA)) then
            call api_error(error%ptr, "get_cavity_gradient", &
                           "Gradient not computed - call compute_cavity_gradient first")
            return
         end if

         nsph = cavity%nsph
         ngrid = cavity%ngrid

         ! Validate the caller's capacities
         if (nsph_cap < nsph .or. ngrid_cap < ngrid) then
            call api_error(error%ptr, "get_cavity_gradient", &
                           "Array capacity is too small - use get_cavity_sizes")
            return
         end if

         ! Copy gradient arrays
         A_tot1_rA(:, :nsph) = cavity%A_tot1_rA(:, :)
         V_tot1_rA(:, :nsph) = cavity%V_tot1_rA(:, :)
         asph1_rA(:, :nsph, :nsph) = cavity%asph1_rA(:, :, :)
         vsph1_rA(:, :nsph, :nsph) = cavity%vsph1_rA(:, :, :)
         xyz1_rA(:, :, :nsph, :ngrid) = cavity%xyz1_rA(:, :, :, :)
         r_iI1_rA(:, :nsph, :ngrid) = cavity%r_iI1_rA(:, :, :)
         rho1_rA(:, :nsph, :ngrid) = cavity%rho1_rA(:, :, :)

      class default
         call api_error(error%ptr, "get_cavity_gradient", &
                        "Cavity is not DROP type - gradient API only supports DROP cavities")
      end select

   end subroutine get_cavity_gradient_api

   !> Assemble the Gaussian PCM A-matrix together with its nuclear derivatives
   !>
   !> - call compute_cavity_gradient first for the gradient
   !> - arrays in Fortran shape notation (column-major: the leftmost index is
   !>   contiguous); a C caller passes flat buffers of the same total size:
   !>   Amat0(ngrid, ngrid)                  - A-matrix (symmetric)
   !>   Amat1_rA(3, nsph, ngrid, ngrid)      - gradient of A-matrix
   !>   xi(ngrid)                            - xi values
   subroutine get_amat_gradient_api(verror, vcav, nsph_cap, ngrid_cap, &
         & Amat0, Amat1_rA, xi) &
         & bind(C, name=namespace//"get_amat_gradient")
      !> Required diagnostic handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Allocated sphere capacity
      integer(c_int), value :: nsph_cap
      !> Allocated grid-point capacity
      integer(c_int), value :: ngrid_cap
      !> PCM interaction matrix
      real(c_double), intent(inout), optional :: Amat0(ngrid_cap, ngrid_cap)
      !> PCM matrix derivatives with respect to nuclear coordinates
      real(c_double), intent(inout), optional :: Amat1_rA(3, nsph_cap, ngrid_cap, ngrid_cap)
      !> Gaussian inverse widths
      real(c_double), intent(inout), optional :: xi(ngrid_cap)
      !> Decoded error handle for argument validation
      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav
      type(error_type), allocatable :: cavity_error
      real(wp), allocatable :: Amat0_f(:, :)
      real(wp), allocatable :: Amat1_rA_f(:, :, :, :)
      integer :: ngrid, nsph

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)
      if (.not. present(Amat0)) then
         call api_error(error%ptr, "get_amat_gradient", "Required pointer 'Amat0' is missing")
         return
      end if
      if (.not. present(Amat1_rA)) then
         call api_error(error%ptr, "get_amat_gradient", "Required pointer 'Amat1_rA' is missing")
         return
      end if
      if (.not. present(xi)) then
         call api_error(error%ptr, "get_amat_gradient", "Required pointer 'xi' is missing")
         return
      end if
      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_amat_gradient", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_amat_gradient", "Cavity is not initialized")
         return
      end if

      associate (cavity => cav%ptr)
         if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f) &
             .or. .not. allocated(cavity%xyz)) then
            call api_error(error%ptr, "get_amat_gradient", &
                           "Cavity does not provide a Gaussian PCM surface")
            return
         end if
         if (.not. allocated(cavity%xi1_rA) .or. .not. allocated(cavity%f1_rA) &
             .or. .not. allocated(cavity%xyz1_rA)) then
            call api_error(error%ptr, "get_amat_gradient", &
                           "Gradient not computed - call compute_cavity_gradient first")
            return
         end if

         ngrid = cavity%ngrid
         nsph = cavity%nsph

         if (ngrid_cap < ngrid .or. nsph_cap < nsph) then
            call api_error(error%ptr, "get_amat_gradient", &
                           "Array capacity is too small - use get_cavity_sizes")
            return
         end if

         allocate (Amat0_f(ngrid, ngrid), Amat1_rA_f(3, nsph, ngrid, ngrid))
         call assemble_pcm_amat_with_gradient(cavity%xi0, cavity%f, cavity%xyz, &
                                              cavity%xi1_rA, cavity%f1_rA, cavity%xyz1_rA, &
                                              Amat0_f, Amat1_rA_f, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "get_amat_gradient", cavity_error%message)
            return
         end if

         Amat0(:ngrid, :ngrid) = Amat0_f(:, :)
         Amat1_rA(:, :nsph, :ngrid, :ngrid) = Amat1_rA_f(:, :, :, :)
         xi(:ngrid) = cavity%xi0(:)
      end associate

   end subroutine get_amat_gradient_api

   !> Contract Gaussian PCM A-matrix derivatives with two grid vectors
   !>
   !> - grad_rA = sum_ij q1_i (dA_ij/dR_A) q2_j
   !> - call compute_cavity_gradient first
   subroutine contract_amat1_q1q2_rA_api(verror, vcav, c_q1, c_q2, c_grad_rA) &
         & bind(C, name=namespace//"contract_amat1_q1q2_rA")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Cavity handle
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      !> First contraction vector, length ngrid
      type(c_ptr), value :: c_q1
      !> Second contraction vector, length ngrid
      type(c_ptr), value :: c_q2
      !> Output nuclear gradient, C shape (natoms,3)
      type(c_ptr), value :: c_grad_rA
      real(c_double), pointer :: q1(:)
      real(c_double), pointer :: q2(:)
      real(c_double), pointer :: grad_rA(:, :)
      type(error_type), allocatable :: cavity_error
      integer :: ngrid, nsph

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "contract_amat1_q1q2_rA", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "contract_amat1_q1q2_rA", "Cavity is not initialized")
         return
      end if

      associate (cavity => cav%ptr)
         if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f) &
             .or. .not. allocated(cavity%xyz) .or. .not. allocated(cavity%xi1_rA) &
             .or. .not. allocated(cavity%f1_rA) .or. .not. allocated(cavity%xyz1_rA)) then
            call api_error(error%ptr, "contract_amat1_q1q2_rA", &
                           "Gaussian PCM derivatives are unavailable - "// &
                           "call compute_cavity_gradient first")
            return
         end if

         ngrid = cavity%ngrid
         nsph = cavity%nsph

         if (.not. c_associated(c_q1) .or. .not. c_associated(c_q2) .or. .not. c_associated(c_grad_rA)) then
            call api_error(error%ptr, "contract_amat1_q1q2_rA", "Null array pointer provided")
            return
         end if

         call c_f_pointer(c_q1, q1, [ngrid])
         call c_f_pointer(c_q2, q2, [ngrid])
         call c_f_pointer(c_grad_rA, grad_rA, [3, nsph])

         block
            !> Gaussian-surface adjoint channels of q1^T dA q2
            real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)

            allocate (w_xi(ngrid), w_f(ngrid), w_xyz(3, ngrid))
            call pcm_amat_surface_weights(cavity%xi0, cavity%f, cavity%xyz, &
                                          q1, q2, w_xi, w_f, w_xyz, cavity_error)
            if (.not. allocated(cavity_error)) then
               call pcm_amat_nuclear_gradient(cavity%xi1_rA, cavity%f1_rA, &
                                              cavity%xyz1_rA, w_xi, w_f, w_xyz, &
                                              grad_rA, cavity_error)
            end if
         end block
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_amat1_q1q2_rA", cavity_error%message)
            return
         end if
      end associate

   end subroutine contract_amat1_q1q2_rA_api

   !> Contract Gaussian PCM A-matrix derivatives to per-grid surface weights
   !>
   !> - weights w_xi, w_f and w_xyz satisfying
   !>   q1^T dA q2 = sum_i w_xi_i dxi_i + w_f_i df_i + w_xyz(:,i).dxyz_i
   subroutine contract_amat1_q1q2_surface_weights_api(verror, vcav, c_q1, c_q2, &
         & c_w_xi, c_w_f, c_w_xyz) &
         & bind(C, name=namespace//"contract_amat1_q1q2_surface_weights")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Cavity handle
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      !> First contraction vector, length ngrid
      type(c_ptr), value :: c_q1
      !> Second contraction vector, length ngrid
      type(c_ptr), value :: c_q2
      !> Gaussian-width weights, length ngrid
      type(c_ptr), value :: c_w_xi
      !> Switching-factor weights, length ngrid
      type(c_ptr), value :: c_w_f
      !> Surface-position weights, C shape (ngrid,3)
      type(c_ptr), value :: c_w_xyz
      real(c_double), pointer :: q1(:)
      real(c_double), pointer :: q2(:)
      real(c_double), pointer :: w_xi(:)
      real(c_double), pointer :: w_f(:)
      real(c_double), pointer :: w_xyz(:, :)
      type(error_type), allocatable :: cavity_error
      integer :: ngrid

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "contract_amat1_q1q2_surface_weights", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "contract_amat1_q1q2_surface_weights", "Cavity is not initialized")
         return
      end if

      associate (cavity => cav%ptr)
         if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f) &
             .or. .not. allocated(cavity%xyz)) then
            call api_error(error%ptr, "contract_amat1_q1q2_surface_weights", &
                           "Cavity does not provide a Gaussian PCM surface")
            return
         end if

         ngrid = cavity%ngrid

         if (.not. c_associated(c_q1) .or. .not. c_associated(c_q2) &
             .or. .not. c_associated(c_w_xi) .or. .not. c_associated(c_w_f) &
             .or. .not. c_associated(c_w_xyz)) then
            call api_error(error%ptr, "contract_amat1_q1q2_surface_weights", "Null array pointer provided")
            return
         end if

         call c_f_pointer(c_q1, q1, [ngrid])
         call c_f_pointer(c_q2, q2, [ngrid])
         call c_f_pointer(c_w_xi, w_xi, [ngrid])
         call c_f_pointer(c_w_f, w_f, [ngrid])
         call c_f_pointer(c_w_xyz, w_xyz, [3, ngrid])

         call pcm_amat_surface_weights(cavity%xi0, cavity%f, cavity%xyz, &
                                       q1, q2, w_xi, w_f, w_xyz, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_amat1_q1q2_surface_weights", cavity_error%message)
            return
         end if
      end associate

   end subroutine contract_amat1_q1q2_surface_weights_api

   !> Contract the original DROP surface channels to LSF adjoint weights
   !>
   !> - preserves the version-0.5 ABI
   !> - use the extended entry point for the normal and principal-curvature
   !>   channels
   subroutine contract_surface_lsf_weights_api(verror, vcav, c_w_xi, c_w_f, c_w_xyz, &
         & c_w_lsf0, c_w_lsf1, c_w_lsf2) &
         & bind(C, name=namespace//"contract_surface_lsf_weights")
      !> Error handle
      type(c_ptr), value :: verror
      !> Cavity handle
      type(c_ptr), value :: vcav
      !> Gaussian-width input weights
      type(c_ptr), value :: c_w_xi
      !> Switching-factor input weights
      type(c_ptr), value :: c_w_f
      !> Surface-position input weights
      type(c_ptr), value :: c_w_xyz
      !> Level-set value output weights
      type(c_ptr), value :: c_w_lsf0
      !> Level-set gradient output weights
      type(c_ptr), value :: c_w_lsf1
      !> Level-set Hessian output weights
      type(c_ptr), value :: c_w_lsf2

      !> Decoded diagnostic handle
      type(vp_error), pointer :: error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      call contract_surface_lsf_weights_extended_api(verror, vcav, c_w_xi, c_w_f, &
                                                     c_w_xyz, c_w_lsf0, c_w_lsf1, c_w_lsf2, &
                                                     c_null_ptr, c_null_ptr, c_null_ptr)
      call prefix_api_error(error%ptr, "contract_surface_lsf_weights")

   end subroutine contract_surface_lsf_weights_api

   !> Contract DROP surface weights to per-grid LSF adjoint weights
   !>
   !> - the projected-coordinate and xi chains are always contracted
   !> - the optional outward-normal (c_w_n) and principal-curvature (c_w_k1,
   !>   c_w_k2) channels are folded in when the caller supplies them; NULL skips
   !>   the channel
   subroutine contract_surface_lsf_weights_extended_api(verror, vcav, c_w_xi, c_w_f, c_w_xyz, &
         & c_w_lsf0, c_w_lsf1, c_w_lsf2, c_w_n, c_w_k1, c_w_k2) &
         & bind(C, name=namespace//"contract_surface_lsf_weights_extended")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Cavity handle
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      !> Gaussian-width weights, length ngrid
      type(c_ptr), value :: c_w_xi
      !> Switching-factor weights, length ngrid
      type(c_ptr), value :: c_w_f
      !> Surface-position weights, C shape (ngrid,3)
      type(c_ptr), value :: c_w_xyz
      !> Output level-set value weights, length ngrid
      type(c_ptr), value :: c_w_lsf0
      !> Output level-set gradient weights, C shape (ngrid,3)
      type(c_ptr), value :: c_w_lsf1
      !> Output level-set Hessian weights, C shape (ngrid,3,3)
      type(c_ptr), value :: c_w_lsf2
      !> Optional surface weights for the outward normal (3, ngrid); NULL to skip
      type(c_ptr), value :: c_w_n
      !> Optional surface weights for the first principal curvature (ngrid); NULL to skip
      type(c_ptr), value :: c_w_k1
      !> Optional surface weights for the second principal curvature (ngrid); NULL to skip
      type(c_ptr), value :: c_w_k2
      real(c_double), pointer :: w_xi(:)
      real(c_double), pointer :: w_f(:)
      real(c_double), pointer :: w_xyz(:, :)
      real(c_double), pointer :: w_lsf0(:)
      real(c_double), pointer :: w_lsf1(:, :)
      real(c_double), pointer :: w_lsf2(:, :, :)
      !> Decoded optional normal weights
      real(c_double), pointer :: w_n(:, :)
      !> Decoded optional curvature weights
      real(c_double), pointer :: w_k1(:), w_k2(:)
      !> Packed surface-adjoint accumulator handed to the cavity
      type(cavity_surface_adjoint_type) :: acc
      type(error_type), allocatable :: cavity_error
      integer :: ngrid

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "contract_surface_lsf_weights_extended", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "contract_surface_lsf_weights_extended", "Cavity is not initialized")
         return
      end if

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         ngrid = cavity%ngrid

         if (.not. c_associated(c_w_xi) .or. .not. c_associated(c_w_f) &
             .or. .not. c_associated(c_w_xyz) .or. .not. c_associated(c_w_lsf0) &
             .or. .not. c_associated(c_w_lsf1) .or. .not. c_associated(c_w_lsf2)) then
            call api_error(error%ptr, "contract_surface_lsf_weights_extended", "Null array pointer provided")
            return
         end if

         call c_f_pointer(c_w_xi, w_xi, [ngrid])
         call c_f_pointer(c_w_f, w_f, [ngrid])
         call c_f_pointer(c_w_xyz, w_xyz, [3, ngrid])
         call c_f_pointer(c_w_lsf0, w_lsf0, [ngrid])
         call c_f_pointer(c_w_lsf1, w_lsf1, [3, ngrid])
         call c_f_pointer(c_w_lsf2, w_lsf2, [3, 3, ngrid])

         ! Optional channels: a NULL pointer leaves the matching accumulator
         ! channel at zero, so the contraction skips that channel entirely
         w_n => null()
         w_k1 => null()
         w_k2 => null()
         if (c_associated(c_w_n)) call c_f_pointer(c_w_n, w_n, [3, ngrid])
         if (c_associated(c_w_k1)) call c_f_pointer(c_w_k1, w_k1, [ngrid])
         if (c_associated(c_w_k2)) call c_f_pointer(c_w_k2, w_k2, [ngrid])

         call acc%init(ngrid)
         call acc%add_surface_weights(cavity_error, w_xi=w_xi, w_f=w_f, w_xyz=w_xyz)
         if (.not. allocated(cavity_error) .and. associated(w_n)) then
            call acc%add_surface_weights(cavity_error, w_n=w_n)
         end if
         if (.not. allocated(cavity_error) .and. associated(w_k1)) then
            call acc%add_surface_weights(cavity_error, w_k1=w_k1)
         end if
         if (.not. allocated(cavity_error) .and. associated(w_k2)) then
            call acc%add_surface_weights(cavity_error, w_k2=w_k2)
         end if
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_surface_lsf_weights_extended", cavity_error%message)
            return
         end if

         call cavity%contract_surface_lsf_weights(acc, w_lsf0, w_lsf1, w_lsf2, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_surface_lsf_weights_extended", cavity_error%message)
            return
         end if
      class default
         call api_error(error%ptr, "contract_surface_lsf_weights_extended", &
                        "Cavity is not DROP type - contraction API only supports DROP cavities")
      end select

   end subroutine contract_surface_lsf_weights_extended_api

   !> Contract the direct nuclear term and the host total surface-position weight
   !>
   !> - the array behind `c_w_phi` is the potential adjoint `dE/dphi_i`, the
   !>   surface charge of a stationary PCM
   !> - the array behind `c_w_xyz` is the host *total* `w_xyz_i = dE_host/dr_i`
   !>   (nuclear plus electronic), contracted with the surface response as is
   !> - cavity-level utility, outside the model coupling protocol
   !> - call compute_cavity_gradient first
   subroutine contract_pcm_nuclear_gradient_api(verror, vcav, c_w_phi, c_w_xyz, c_za, c_grad_rA) &
         & bind(C, name=namespace//"contract_pcm_nuclear_gradient")
      !> Diagnostic handle
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      !> Cavity handle
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      !> Potential adjoint `w_phi = dE/dphi`, length ngrid
      type(c_ptr), value :: c_w_phi
      !> Surface-position weights, C shape (ngrid,3)
      type(c_ptr), value :: c_w_xyz
      !> Nuclear charges, length natoms
      type(c_ptr), value :: c_za
      !> Output nuclear gradient, C shape (natoms,3)
      type(c_ptr), value :: c_grad_rA
      real(c_double), pointer :: w_phi(:)
      real(c_double), pointer :: w_xyz(:, :)
      real(c_double), pointer :: za(:)
      real(c_double), pointer :: grad_rA(:, :)
      type(error_type), allocatable :: cavity_error
      integer :: ngrid, nsph

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)
      if (allocated(error%ptr)) deallocate (error%ptr)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "contract_pcm_nuclear_gradient", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "contract_pcm_nuclear_gradient", "Cavity is not initialized")
         return
      end if

      associate (cavity => cav%ptr)
         if (.not. allocated(cavity%xyz) .or. .not. allocated(cavity%sphxyz) &
             .or. .not. allocated(cavity%xyz1_rA)) then
            call api_error(error%ptr, "contract_pcm_nuclear_gradient", &
                           "Cavity position derivatives are unavailable - "// &
                           "call compute_cavity_gradient first")
            return
         end if

         ngrid = cavity%ngrid
         nsph = cavity%nsph

         if (.not. c_associated(c_w_phi) .or. .not. c_associated(c_w_xyz) &
             .or. .not. c_associated(c_za) .or. .not. c_associated(c_grad_rA)) then
            call api_error(error%ptr, "contract_pcm_nuclear_gradient", "Null array pointer provided")
            return
         end if

         call c_f_pointer(c_w_phi, w_phi, [ngrid])
         call c_f_pointer(c_w_xyz, w_xyz, [3, ngrid])
         call c_f_pointer(c_za, za, [nsph])
         call c_f_pointer(c_grad_rA, grad_rA, [3, nsph])

         call pcm_electrostatic_nuclear_gradient(cavity%xyz, cavity%sphxyz, &
                                                 cavity%xyz1_rA, w_phi, w_xyz, za, &
                                                 grad_rA, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_pcm_nuclear_gradient", cavity_error%message)
            return
         end if
      end associate

   end subroutine contract_pcm_nuclear_gradient_api

   !> Delete a cavity of any concrete type
   subroutine delete_cavity_api(vcav) &
         & bind(C, name=namespace//"delete_cavity")
      !> Cavity handle
      type(c_ptr), intent(inout), optional :: vcav
      type(vp_cavity), pointer :: cav

      if (.not. present(vcav)) return
      if (c_associated(vcav)) then
         call c_f_pointer(vcav, cav)
         if (cav%owned .and. associated(cav%ptr)) deallocate (cav%ptr)
         nullify (cav%ptr)
         call cav%ctx%delete()
         deallocate (cav)
         vcav = c_null_ptr
      end if

   end subroutine delete_cavity_api

   !> Copy a Fortran string into a bounded, null-terminated C buffer
   subroutine f_c_character(rhs, lhs, len)
      !> Destination C buffer, including the terminator
      character(kind=c_char), intent(out) :: lhs(*)
      !> Source Fortran string
      character(len=*), intent(in) :: rhs
      !> Allocated destination capacity in characters
      integer, intent(in) :: len
      integer :: length

      if (len <= 0) return

      length = min(len - 1, len_trim(rhs))

      if (length > 0) lhs(1:length) = transfer(rhs(1:length), lhs(1:length))
      lhs(length + 1:length + 1) = c_null_char

   end subroutine f_c_character

   !> Decode a bounded C string; report whether a terminator was found
   subroutine c_f_character_ptr(rhs_ptr, lhs, max_len, truncated)
      !> Optional source C string pointer
      type(c_ptr), value, intent(in) :: rhs_ptr
      !> Decoded string; empty for a NULL source
      character(len=:, kind=c_char), allocatable, intent(out) :: lhs
      !> Maximum number of characters to inspect
      integer, intent(in) :: max_len
      !> Whether the source lacks a terminator within max_len
      logical, intent(out), optional :: truncated
      character(kind=c_char), pointer :: rhs(:)
      integer :: ii, nchar, scan_len
      logical :: has_null

      if (present(truncated)) truncated = .false.

      if (.not. c_associated(rhs_ptr)) then
         !> Decoded string; empty for a NULL source
         allocate (character(len=0, kind=c_char) :: lhs)
         return
      end if

      scan_len = max(1, min(max_len, huge(scan_len) - 1))
      call c_f_pointer(rhs_ptr, rhs, [scan_len])

      has_null = .false.
      do ii = 1, scan_len
         if (rhs(ii) == c_null_char) then
            has_null = .true.
            exit
         end if
      end do

      if (has_null) then
         nchar = ii - 1
      else
         nchar = scan_len
         if (present(truncated)) truncated = .true.
      end if

      !> Decoded string; empty for a NULL source
      allocate (character(len=nchar, kind=c_char) :: lhs)
      if (nchar > 0) lhs = transfer(rhs(1:nchar), lhs)

   end subroutine c_f_character_ptr

   !> Reject coincident atoms before accepting a molecular geometry
   subroutine verify_structure(error, mol)
      !> Diagnostic on failure
      type(error_type), allocatable, intent(out) :: error
      !> Molecular geometry in Bohr
      type(structure_type), intent(in) :: mol
      integer :: iat, jat, stat
      stat = 0
      do iat = 1, mol%nat
         do jat = 1, iat - 1
            if (norm2(mol%xyz(:, jat) - mol%xyz(:, iat)) < 1.0e-9_wp) stat = stat + 1
         end do
      end do
      if (stat > 0) then
         call fatal_error(error, "Too close interatomic distances found")
      end if
   end subroutine verify_structure

end module moist_api
