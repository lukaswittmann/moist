!> Defensive-preamble tests for the C API entry points
module test_api
   use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_null_ptr, c_int, c_double, c_bool, &
                                          c_funptr, c_funloc, c_associated, c_f_pointer, &
                                          c_char, c_null_char, c_size_t, c_sizeof, c_null_funptr, &
                                          c_int8_t
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io_structure, only: structure_type
   use moist_api, only: vp_cavity, vp_error, vp_response, response_get_api, &
      & next_response_item_api, response_item_name_api, vp_structure, vp_radii, &
      & vp_model, update_structure_api, set_custom_radii_atoms_api, &
      & set_custom_radii_elements_api, update_solvation_model_api, &
      & get_solvation_model_cavity_api, delete_solvation_model_api, &
      & new_pv_component_api, delete_solvation_component_api, &
      & general_model_add_component_api, new_coupling_api, delete_coupling_api, &
      & new_response_api, delete_response_api, general_model_get_gradient_api, &
      & set_isodensity_density_api, get_drop_cavity_tolerance_api, &
      & compute_cavity_gradient_api, compute_anchor_gradient_api, &
      & get_anchor_gradient_api, get_cavity_gradient_api, get_amat_gradient_api, &
      & contract_amat1_q1q2_surface_weights_api, &
      & contract_surface_lsf_weights_extended_api, contract_pcm_nuclear_gradient_api
   use moist_cavity_fields, only: cavity_field_query_type
   use moist_cavity_type, only: cavity_type
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: density_response_type, &
      & potential_adjoint_response_type, response_accumulate, response_type
   use moist_model_type, only: solvation_model_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_api

   !> Frozen isodensity options layout exercised through the C interface
   type, bind(C) :: iso_options
      !> Caller layout size
      integer(c_size_t) :: struct_size
      !> Density isovalue
      real(c_double) :: rho_iso
      !> Level-set scale
      real(c_double) :: scale
   end type iso_options

   !> Frozen DROP options layout exercised through the C interface
   type, bind(C) :: drop_options
      !> Caller layout size
      integer(c_size_t) :: struct_size
      !> Lebedev points per atom
      integer(c_int) :: nleb
      !> Diagnostic checks
      logical(c_bool) :: debug
      !> Output detail level
      integer(c_int) :: verbosity
      !> Fine cavity refinement
      logical(c_bool) :: do_fine
      !> Projection convergence tolerance
      real(c_double) :: tolerance
      !> Maximum projection iterations
      integer(c_int) :: proj_maxiter
      !> Projection refinement level
      integer(c_int) :: proj_level
      !> Branch-selection smoothing weight
      real(c_double) :: branch_weight_s
      !> Grid density kernel length
      real(c_double) :: rho_grid_h
      !> Lebedev pruning level
      integer(c_int) :: wleb_prune_level
      !> Reserved ABI storage
      integer(c_int) :: reserved0
   end type drop_options

   !> Test double: a cavity whose arrays the test fills by hand
   !>
   !> Reaches validation branches no built-in cavity can: a width that is not
   !> positive, extents that disagree with each other, and a field name longer
   !> than the C name buffer
   type, extends(cavity_type) :: stub_cavity
      !> Name of the one field the stub declares
      character(len=:), allocatable :: field_name
   contains
      !> Accept any structure without building anything
      procedure :: update => stub_cavity_update
      !> Report no derivatives
      procedure :: get_gradient => stub_cavity_gradient
      !> Declare `field_name` over the widths
      procedure :: list_fields => stub_cavity_fields
   end type stub_cavity

   !> Test double: a solvation model that is not the general model
   type, extends(solvation_model_type) :: stub_model
   contains
      !> Accept any structure
      procedure :: update => stub_model_update
      !> Contribute no energy
      procedure :: get_energy => stub_model_energy
      !> Contribute no response
      procedure :: get_response => stub_model_response
      !> Contribute no gradient
      procedure :: get_gradient => stub_model_gradient
   end type stub_model

   !> C bindings of the entry points under test
   interface

      !> C binding for moist_new_svdw_lsf
      !>
      !> @param[in] error Native error
      !> @param[in] options Native options
      function moist_new_svdw_lsf(error, options) result(lsf) bind(C)
         import :: c_ptr
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: error
         !> Native options
         type(c_ptr), value, intent(in) :: options
         !> Owning result handle
         type(c_ptr) :: lsf
      end function moist_new_svdw_lsf

      !> C binding for moist_new_isodensity_callback_lsf
      !>
      !> @param[in] error Native error
      !> @param[in] callback Native callback
      !> @param[in] context Native context
      !> @param[in] options Native options
      function moist_new_isodensity_callback_lsf(error, callback, context, options) result(lsf) bind(C)
         import :: c_ptr, c_funptr
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: error
         !> Native callback
         type(c_funptr), value, intent(in) :: callback
         !> Native context
         type(c_ptr), value, intent(in) :: context
         !> Native options
         type(c_ptr), value, intent(in) :: options
         !> Owning result handle
         type(c_ptr) :: lsf
      end function moist_new_isodensity_callback_lsf

      !> C binding for moist_new_drop_cavity
      !>
      !> @param[in] error Native error
      !> @param[in] lsf Native lsf
      !> @param[in] radii Native radii
      !> @param[in] options Native options
      function moist_new_drop_cavity(error, lsf, radii, options) result(cavity) bind(C)
         import :: c_ptr
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: error
         !> Native lsf
         type(c_ptr), value, intent(in) :: lsf
         !> Native radii
         type(c_ptr), value, intent(in) :: radii
         !> Native options
         type(c_ptr), value, intent(in) :: options
         !> Owning result handle
         type(c_ptr) :: cavity
      end function moist_new_drop_cavity

      !> C binding for moist_delete_lsf
      !>
      !> @param[in,out] lsf Native lsf
      subroutine moist_delete_lsf(lsf) bind(C)
         import :: c_ptr
         implicit none(type, external)
         !> Native lsf
         type(c_ptr), intent(inout) :: lsf
      end subroutine moist_delete_lsf

      subroutine moist_get_cavity_sizes(verror, vcav, ngrid, nsph) &
            & bind(C, name="moist_get_cavity_sizes")
         import :: c_ptr, c_int
         implicit none(type, external)
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), intent(inout), optional :: ngrid
         integer(c_int), intent(inout), optional :: nsph
      end subroutine moist_get_cavity_sizes

      function moist_new_structure(verror, natoms, numbers, positions, &
            & lattice, periodic) result(vmol) &
            & bind(C, name="moist_new_structure")
         import :: c_ptr, c_int, c_double, c_bool
         implicit none(type, external)
         type(c_ptr), value :: verror
         integer(c_int), value :: natoms
         integer(c_int), intent(in), optional :: numbers(natoms)
         real(c_double), intent(in), optional :: positions(3, natoms)
         real(c_double), intent(in), optional :: lattice(3, 3)
         logical(c_bool), intent(in), optional :: periodic(3)
         type(c_ptr) :: vmol
      end function moist_new_structure

      subroutine moist_delete_structure(vmol) &
            & bind(C, name="moist_delete_structure")
         import :: c_ptr
         implicit none(type, external)
         type(c_ptr), intent(inout) :: vmol
      end subroutine moist_delete_structure


      subroutine moist_update_cavity(verror, vcav, vmol) &
            & bind(C, name="moist_update_cavity")
         import :: c_ptr
         implicit none(type, external)
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         type(c_ptr), value :: vmol
      end subroutine moist_update_cavity

      subroutine moist_delete_cavity(vcav) &
            & bind(C, name="moist_delete_cavity")
         import :: c_ptr
         implicit none(type, external)
         type(c_ptr), intent(inout) :: vcav
      end subroutine moist_delete_cavity

      subroutine moist_get_cavity_results(verror, vcav, ngrid_cap, nsph_cap, &
            & area, volume, ngrid, nsph, xyz, a, owner, converged, radii, asph) &
            & bind(C, name="moist_get_cavity_results")
         import :: c_ptr, c_int, c_double, c_bool
         implicit none(type, external)
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ngrid_cap
         integer(c_int), value :: nsph_cap
         real(c_double), intent(inout), optional :: area
         real(c_double), intent(inout), optional :: volume
         integer(c_int), intent(inout), optional :: ngrid
         integer(c_int), intent(inout), optional :: nsph
         real(c_double), intent(inout), optional :: xyz(3, ngrid_cap)
         real(c_double), intent(inout), optional :: a(ngrid_cap)
         integer(c_int), intent(inout), optional :: owner(ngrid_cap)
         logical(c_bool), intent(inout), optional :: converged(ngrid_cap)
         real(c_double), intent(inout), optional :: radii(nsph_cap)
         real(c_double), intent(inout), optional :: asph(nsph_cap)
      end subroutine moist_get_cavity_results

      subroutine moist_get_cavity_field_count(verror, vcav, nfield) &
            & bind(C, name="moist_get_cavity_field_count")
         import :: c_ptr, c_int
         implicit none(type, external)
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), intent(inout), optional :: nfield
      end subroutine moist_get_cavity_field_count

      subroutine moist_get_cavity_field_info(verror, vcav, ifield, name, &
            & dtype, rank, dims, count) &
            & bind(C, name="moist_get_cavity_field_info")
         import :: c_ptr, c_int, c_char
         implicit none(type, external)
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ifield
         character(kind=c_char), intent(inout), optional :: name(*)
         integer(c_int), intent(inout), optional :: dtype
         integer(c_int), intent(inout), optional :: rank
         integer(c_int), intent(inout), optional :: dims(2)
         integer(c_int), intent(inout), optional :: count
      end subroutine moist_get_cavity_field_info

      subroutine moist_get_cavity_field_real(verror, vcav, cname, values) &
            & bind(C, name="moist_get_cavity_field_real")
         import :: c_ptr, c_double
         implicit none(type, external)
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         type(c_ptr), value :: cname
         real(c_double), intent(inout), optional :: values(*)
      end subroutine moist_get_cavity_field_real

      subroutine moist_get_cavity_field_int(verror, vcav, cname, values) &
            & bind(C, name="moist_get_cavity_field_int")
         import :: c_ptr, c_int
         implicit none(type, external)
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         type(c_ptr), value :: cname
         integer(c_int), intent(inout), optional :: values(*)
      end subroutine moist_get_cavity_field_int

      subroutine moist_get_cavity_field_bool(verror, vcav, cname, values) &
            & bind(C, name="moist_get_cavity_field_bool")
         import :: c_ptr, c_bool
         implicit none(type, external)
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         type(c_ptr), value :: cname
         logical(c_bool), intent(inout), optional :: values(*)
      end subroutine moist_get_cavity_field_bool

      subroutine moist_get_cavity_field_about(verror, vcav, cname, about, capacity, length) &
            & bind(C, name="moist_get_cavity_field_about")
         import :: c_ptr, c_size_t, c_char
         implicit none(type, external)
         type(c_ptr), value, intent(in) :: verror
         type(c_ptr), value, intent(in) :: vcav
         type(c_ptr), value, intent(in) :: cname
         character(kind=c_char), intent(inout), optional :: about(*)
         integer(c_size_t), value, intent(in) :: capacity
         integer(c_size_t), intent(inout), optional :: length
      end subroutine moist_get_cavity_field_about

      subroutine moist_assemble_amat(verror, vcav, ngrid_cap, amat0, xi) &
            & bind(C, name="moist_assemble_amat")
         import :: c_ptr, c_int, c_double
         implicit none(type, external)
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ngrid_cap
         real(c_double), intent(inout), optional :: amat0(ngrid_cap, ngrid_cap)
         real(c_double), intent(inout), optional :: xi(ngrid_cap)
      end subroutine moist_assemble_amat

      !> C binding for moist_init_drop_options
      !>
      !> @param[in] verror Native error
      !> @param[in] options Native options buffer
      !> @param[in] bytes Native buffer size
      subroutine moist_init_drop_options(verror, options, bytes) bind(C)
         import :: c_ptr, c_size_t
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native options buffer
         type(c_ptr), value, intent(in) :: options
         !> Native buffer size
         integer(c_size_t), value, intent(in) :: bytes
      end subroutine moist_init_drop_options

      !> C binding for moist_init_iswig_options
      !>
      !> @param[in] verror Native error
      !> @param[in] options Native options buffer
      !> @param[in] bytes Native buffer size
      subroutine moist_init_iswig_options(verror, options, bytes) bind(C)
         import :: c_ptr, c_size_t
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native options buffer
         type(c_ptr), value, intent(in) :: options
         !> Native buffer size
         integer(c_size_t), value, intent(in) :: bytes
      end subroutine moist_init_iswig_options

      !> C binding for moist_init_svdw_options
      !>
      !> @param[in] verror Native error
      !> @param[in] options Native options buffer
      !> @param[in] bytes Native buffer size
      subroutine moist_init_svdw_options(verror, options, bytes) bind(C)
         import :: c_ptr, c_size_t
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native options buffer
         type(c_ptr), value, intent(in) :: options
         !> Native buffer size
         integer(c_size_t), value, intent(in) :: bytes
      end subroutine moist_init_svdw_options

      !> C binding for moist_init_cfc_options
      !>
      !> @param[in] verror Native error
      !> @param[in] options Native options buffer
      !> @param[in] bytes Native buffer size
      subroutine moist_init_cfc_options(verror, options, bytes) bind(C)
         import :: c_ptr, c_size_t
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native options buffer
         type(c_ptr), value, intent(in) :: options
         !> Native buffer size
         integer(c_size_t), value, intent(in) :: bytes
      end subroutine moist_init_cfc_options

      !> C binding for moist_init_isodensity_options
      !>
      !> @param[in] verror Native error
      !> @param[in] options Native options buffer
      !> @param[in] bytes Native buffer size
      subroutine moist_init_isodensity_options(verror, options, bytes) bind(C)
         import :: c_ptr, c_size_t
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native options buffer
         type(c_ptr), value, intent(in) :: options
         !> Native buffer size
         integer(c_size_t), value, intent(in) :: bytes
      end subroutine moist_init_isodensity_options

      !> C binding for moist_init_model_options
      !>
      !> @param[in] verror Native error
      !> @param[in] options Native options buffer
      !> @param[in] bytes Native buffer size
      subroutine moist_init_model_options(verror, options, bytes) bind(C)
         import :: c_ptr, c_size_t
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native options buffer
         type(c_ptr), value, intent(in) :: options
         !> Native buffer size
         integer(c_size_t), value, intent(in) :: bytes
      end subroutine moist_init_model_options

      !> C binding for moist_init_pcm_options
      !>
      !> @param[in] verror Native error
      !> @param[in] options Native options buffer
      !> @param[in] bytes Native buffer size
      subroutine moist_init_pcm_options(verror, options, bytes) bind(C)
         import :: c_ptr, c_size_t
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native options buffer
         type(c_ptr), value, intent(in) :: options
         !> Native buffer size
         integer(c_size_t), value, intent(in) :: bytes
      end subroutine moist_init_pcm_options

      !> C binding for moist_new_isodensity_lsf
      !>
      !> @param[in] verror Native error
      !> @param[in] nshell Native shell count
      !> @param[in] shell_atom Native shell atoms
      !> @param[in] shell_l Native shell angular momenta
      !> @param[in] shell_nprim Native primitive counts
      !> @param[in] exps Native exponents
      !> @param[in] coeffs Native coefficients
      !> @param[in] options Native options
      function moist_new_isodensity_lsf(verror, nshell, shell_atom, shell_l, shell_nprim, &
            & exps, coeffs, options) result(lsf) bind(C)
         import :: c_ptr, c_int
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native shell count
         integer(c_int), value, intent(in) :: nshell
         !> Native shell atoms
         type(c_ptr), value, intent(in) :: shell_atom
         !> Native shell angular momenta
         type(c_ptr), value, intent(in) :: shell_l
         !> Native primitive counts
         type(c_ptr), value, intent(in) :: shell_nprim
         !> Native exponents
         type(c_ptr), value, intent(in) :: exps
         !> Native coefficients
         type(c_ptr), value, intent(in) :: coeffs
         !> Native options
         type(c_ptr), value, intent(in) :: options
         !> Owning result handle
         type(c_ptr) :: lsf
      end function moist_new_isodensity_lsf

      !> C binding for moist_new_iswig_cavity
      !>
      !> @param[in] verror Native error
      !> @param[in] radii Native radii
      !> @param[in] options Native options
      function moist_new_iswig_cavity(verror, radii, options) result(cavity) bind(C)
         import :: c_ptr
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native radii
         type(c_ptr), value, intent(in) :: radii
         !> Native options
         type(c_ptr), value, intent(in) :: options
         !> Owning result handle
         type(c_ptr) :: cavity
      end function moist_new_iswig_cavity

      !> C binding for moist_new_model
      !>
      !> @param[in] verror Native error
      !> @param[in] cavity Native cavity
      !> @param[in] options Native options
      function moist_new_model(verror, cavity, options) result(model) bind(C)
         import :: c_ptr
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native cavity
         type(c_ptr), value, intent(in) :: cavity
         !> Native options
         type(c_ptr), value, intent(in) :: options
         !> Owning result handle
         type(c_ptr) :: model
      end function moist_new_model

      !> C binding for moist_set_model_isodensity_density
      !>
      !> @param[in] verror Native error
      !> @param[in] model Native model
      !> @param[in] ncart Native density dimension
      !> @param[in] density Native density buffer
      subroutine moist_set_model_isodensity_density(verror, model, ncart, density) bind(C)
         import :: c_ptr, c_int
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: verror
         !> Native model
         type(c_ptr), value, intent(in) :: model
         !> Native density dimension
         integer(c_int), value, intent(in) :: ncart
         !> Native density buffer
         type(c_ptr), value, intent(in) :: density
      end subroutine moist_set_model_isodensity_density

   end interface

   !* ================================================================================= *!
   !*                    Isodensity callback state (failure-channel tests)              *!
   !* ================================================================================= *!

   !> Water geometry the test callback's model density is centered on (Bohr)
   real(c_double), parameter :: cb_centers(3, 3) = reshape( &
      [0.0_c_double, 0.0_c_double, 0.1173_c_double, &
       0.0_c_double, 1.4309_c_double, -0.9370_c_double, &
       0.0_c_double, -1.4309_c_double, -0.9370_c_double], [3, 3])
   !> Exponent of the per-atom s-Gaussian density, 2*alpha for alpha = 0.3
   real(c_double), parameter :: cb_a = 0.6_c_double
   !> Prefactor of the per-atom s-Gaussian density, 2*((2*alpha/pi)**0.75)**2
   real(c_double), parameter :: cb_c = 0.166929_c_double
   !> Isovalue defining the surface (Bohr^-3); puts it ~2.9 Bohr off each atom
   real(c_double), parameter :: cb_rho_iso = 1.0e-3_c_double

   !> Per-test callback state
   !>
   !> Lives in the callback's own `context` rather than in module variables:
   !> test-drive runs the tests of a suite in an OpenMP parallel loop, so two
   !> tests using the same callback would otherwise share one counter
   type :: iso_cb_ctx
      !> Evaluations answered before the callback starts reporting failure
      integer :: fail_after = huge(1)
      !> Status the callback reports once it starts failing
      integer(c_int) :: fail_status = 0_c_int
      !> Evaluation counter, bumped atomically: moist calls the callback from
      !> several OpenMP threads within one build
      integer :: calls = 0
   end type iso_cb_ctx

contains

   !> Collect all C API guard tests
   subroutine collect_api(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("response_arrays", test_response_arrays), &
                  new_unittest("response_walk", test_response_walk), &
                  new_unittest("update_drop_cavity_null_handle", test_update_drop_null_handle), &
                  new_unittest("update_drop_cavity_null_inner", test_update_drop_null_inner), &
                  new_unittest("get_drop_results_null_handle", test_drop_results_null_handle), &
                  new_unittest("get_drop_results_null_inner", test_drop_results_null_inner), &
                  new_unittest("get_cavity_sizes_null_handle", test_cavity_sizes_null_handle), &
                  new_unittest("get_cavity_sizes_null_inner", test_cavity_sizes_null_inner), &
                  new_unittest("cavity_field_null_name", test_cavity_field_null_name), &
                  new_unittest("cavity_field_null_handle", test_cavity_field_null_handle), &
                  new_unittest("cavity_field_null_inner", test_cavity_field_null_inner), &
                  new_unittest("cavity_field_empty_name", test_cavity_field_empty_name), &
                  new_unittest("cavity_field_index_out_of_range", test_cavity_field_index_range), &
                  new_unittest("cavity_field_about_length_and_truncation", test_cavity_field_about_buffer), &
                  new_unittest("cavity_field_bool_accessor", test_cavity_field_bool), &
                  new_unittest("array_capacity_too_small", test_capacity_too_small), &
                  new_unittest("array_capacity_oversized", test_capacity_oversized), &
                  new_unittest("isodensity_callback_fails_first_call", test_iso_callback_fails_first), &
                  new_unittest("isodensity_callback_fails_mid_loop", test_iso_callback_fails_mid_loop), &
                  new_unittest("init_options_without_error_handle", test_init_options_null_error), &
                  new_unittest("lsf_constructor_guards", test_lsf_constructor_guards), &
                  new_unittest("cavity_constructor_guards", test_cavity_constructor_guards), &
                  new_unittest("structure_and_radii_guards", test_structure_radii_guards), &
                  new_unittest("missing_pointers_results", test_missing_pointers_results), &
                  new_unittest("missing_pointers_gradients", test_missing_pointers_gradients), &
                  new_unittest("unbuilt_cavity_guards", test_unbuilt_cavity_guards), &
                  new_unittest("drop_only_entry_points", test_drop_only_entry_points), &
                  new_unittest("contraction_guards", test_contraction_guards), &
                  new_unittest("stub_cavity_guards", test_stub_cavity_guards), &
                  new_unittest("model_handle_guards", test_model_handle_guards), &
                  new_unittest("model_phase_guards", test_model_phase_guards), &
                  new_unittest("response_unavailable_arrays", test_response_unavailable_arrays), &
                  new_unittest("gradient_capacity_too_small", test_gradient_capacity_too_small), &
                  new_unittest("amat_surface_weights", test_amat_surface_weights) &
                  ]

   end subroutine collect_api

   !> Assert that the API error handle carries a message containing `needle`
   subroutine check_api_error(error, err, needle)
      !> Test-drive error
      type(error_type), allocatable, intent(out) :: error
      !> API error handle produced by the entry point
      type(vp_error), pointer, intent(in) :: err
      !> Expected message fragment
      character(len=*), intent(in) :: needle

      if (.not. allocated(err%ptr)) then
         call test_failed(error, "API entry point did not report an error")
         return
      end if
      if (index(err%ptr%message, needle) <= 0) then
         call test_failed(error, "Unexpected API error message: "//err%ptr%message)
         return
      end if

   end subroutine check_api_error

   !> A missing DROP cavity handle is reported, not dereferenced
   subroutine test_update_drop_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err

      allocate (err)
      call moist_update_cavity(c_loc(err), c_null_ptr, c_null_ptr)
      call check_api_error(error, err, "Cavity handle is missing")
      deallocate (err)

   end subroutine test_update_drop_null_handle

   !> A DROP cavity handle with a null inner pointer is reported, not dereferenced
   subroutine test_update_drop_null_inner(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav

      allocate (err)
      allocate (cav)
      call moist_update_cavity(c_loc(err), c_loc(cav), c_null_ptr)
      call check_api_error(error, err, "Cavity is not initialized")
      deallocate (cav)
      deallocate (err)

   end subroutine test_update_drop_null_inner

   !> A missing DROP cavity handle is reported, not dereferenced
   subroutine test_drop_results_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err

      allocate (err)
      call call_get_drop_results(c_loc(err), c_null_ptr)
      call check_api_error(error, err, "Cavity handle is missing")
      deallocate (err)

   end subroutine test_drop_results_null_handle

   !> A DROP cavity handle with a null inner pointer is reported, not dereferenced
   subroutine test_drop_results_null_inner(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav

      allocate (err)
      allocate (cav)
      call call_get_drop_results(c_loc(err), c_loc(cav))
      call check_api_error(error, err, "Cavity is not initialized")
      deallocate (cav)
      deallocate (err)

   end subroutine test_drop_results_null_inner

   !> A missing cavity handle is reported, not dereferenced
   subroutine test_cavity_sizes_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      integer(c_int) :: ngrid, nsph

      allocate (err)
      call moist_get_cavity_sizes(c_loc(err), c_null_ptr, ngrid, nsph)
      call check_api_error(error, err, "Cavity handle is missing")
      deallocate (err)

   end subroutine test_cavity_sizes_null_handle

   !> A cavity handle with a null inner pointer is reported, not dereferenced
   subroutine test_cavity_sizes_null_inner(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav
      integer(c_int) :: ngrid, nsph

      allocate (err)
      allocate (cav)
      call moist_get_cavity_sizes(c_loc(err), c_loc(cav), ngrid, nsph)
      call check_api_error(error, err, "Cavity is not initialized")
      deallocate (cav)
      deallocate (err)

   end subroutine test_cavity_sizes_null_inner

   !> Provide the caller-owned result buffers moist_get_cavity_results writes into
   !> The guarded entry point returns before touching them; they only have to
   !> exist so the call is well formed
   subroutine call_get_drop_results(verror, vcav)
      !> Opaque error handle
      type(c_ptr), intent(in) :: verror
      !> Opaque cavity handle under test
      type(c_ptr), intent(in) :: vcav
      !> Scalar result buffers
      real(c_double) :: area, volume
      integer(c_int) :: ngrid, nsph
      !> Per-point result buffers
      real(c_double) :: xyz(3, 1)
      real(c_double) :: a(1)
      integer(c_int) :: owner(1)
      logical(c_bool) :: converged(1)
      !> Per-sphere result buffers
      real(c_double) :: radii(1), asph(1)

      call moist_get_cavity_results(verror, vcav, 1_c_int, 1_c_int, &
                                    area, volume, ngrid, nsph, xyz, a, owner, converged, radii, asph)

   end subroutine call_get_drop_results

   !> Build a DROP cavity for water through the public C entry points, exactly
   !> as a host would.  The caller owns the returned handles
   subroutine build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      !> Test-drive error
      type(error_type), allocatable, intent(out) :: error
      !> API error handle (allocated here, deallocated by the caller)
      type(vp_error), pointer, intent(out) :: err
      !> Opaque view of `err`
      type(c_ptr), intent(out) :: verror
      !> Structure and cavity handles
      type(c_ptr), intent(out) :: vmol, vcav
      !> Sizes of the built cavity
      integer(c_int), intent(out) :: ngrid, nsph
      !> Water, coordinates in Bohr
      integer(c_int), parameter :: numbers(3) = [8, 1, 1]
      real(c_double), parameter :: positions(3, 3) = reshape( &
         [0.0_c_double, 0.0_c_double, 0.1173_c_double, &
          0.0_c_double, 1.4309_c_double, -0.9370_c_double, &
          0.0_c_double, -1.4309_c_double, -0.9370_c_double], [3, 3])
      !> Lebedev order, kept small so the build stays cheap
      integer(c_int), target :: nleb

      ngrid = 0
      nsph = 0
      nleb = 26_c_int

      allocate (err)
      verror = c_loc(err)

      vmol = moist_new_structure(verror, 3_c_int, numbers, positions)
      block
         type(c_ptr) :: lsf
         lsf = moist_new_svdw_lsf(verror, c_null_ptr)
         vcav = moist_new_drop_cavity(verror, lsf, c_null_ptr, c_null_ptr)
         call moist_delete_lsf(lsf)
      end block
      call moist_update_cavity(verror, vcav, vmol)
      call moist_get_cavity_sizes(verror, vcav, ngrid, nsph)

      if (allocated(err%ptr)) then
         call test_failed(error, "Cavity setup failed: "//err%ptr%message)
         return
      end if
      if (ngrid <= 1 .or. nsph /= 3) then
         call test_failed(error, "Unexpected cavity sizes for water")
         return
      end if

   end subroutine build_water_cavity

   !> Release the handles built by build_water_cavity
   subroutine drop_water_cavity(err, vmol, vcav)
      type(vp_error), pointer, intent(inout) :: err
      type(c_ptr), intent(inout) :: vmol, vcav

      call moist_delete_cavity(vcav)
      call moist_delete_structure(vmol)
      deallocate (err)

   end subroutine drop_water_cavity

   !> Null-terminate a Fortran string into a buffer c_loc can be taken of
   function c_string(text) result(buffer)
      !> Text to hand to a C entry point
      character(len=*), intent(in) :: text
      !> Null-terminated copy
      character(kind=c_char), allocatable :: buffer(:)
      integer :: ichar

      allocate (buffer(len(text) + 1))
      do ichar = 1, len(text)
         buffer(ichar) = text(ichar:ichar)
      end do
      buffer(len(text) + 1) = c_null_char

   end function c_string

   !> A null field name is reported rather than dereferenced
   !>
   !> The name crosses as a bare `const char*` that the entry point has to scan,
   !> so a caller that passes NULL must get an error instead of a walk off
   !> address zero
   subroutine test_cavity_field_null_name(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph
      real(c_double) :: values(1)

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      values = -12345.0_c_double
      call moist_get_cavity_field_real(verror, vcav, c_null_ptr, values)
      call check_api_error(error, err, "Field name is missing")

      if (.not. allocated(error)) then
         if (any(values /= -12345.0_c_double)) then
            call test_failed(error, "get_cavity_field_real wrote for a null name")
         end if
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_cavity_field_null_name

   !* ================================================================================= *!
   !*                          Named field API defensive guards                         *!
   !* ================================================================================= *!

   !> Position of the null terminator in a C string buffer, i.e. its length
   pure integer function c_string_len(buffer) result(nchar)
      !> Null-terminated buffer written by a C entry point
      character(kind=c_char), intent(in) :: buffer(:)
      integer :: ichar

      nchar = size(buffer)
      do ichar = 1, size(buffer)
         if (buffer(ichar) == c_null_char) then
            nchar = ichar - 1
            return
         end if
      end do

   end function c_string_len

   !> A missing cavity handle is reported by the field API, not dereferenced
   !>
   !> `resolve_field_cavity` is the shared preamble of all six field entry
   !> points, so this single guard is what stands between every named read and a
   !> walk off a null handle.  Both an enumerating and a fetching entry point
   !> are exercised because they reach the guard by different routes
   subroutine test_cavity_field_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      character(kind=c_char), allocatable, target :: name_wleb(:)
      integer(c_int) :: nfield
      real(c_double) :: values(1)

      allocate (err)
      name_wleb = c_string("wleb")

      nfield = -1_c_int
      call moist_get_cavity_field_count(c_loc(err), c_null_ptr, nfield)
      call check_api_error(error, err, "Cavity handle is missing")
      if (.not. allocated(error)) then
         if (nfield /= -1_c_int) then
            call test_failed(error, "get_cavity_field_count changed nfield on failure")
         end if
      end if
      if (allocated(err%ptr)) deallocate (err%ptr)

      if (.not. allocated(error)) then
         values = -12345.0_c_double
         call moist_get_cavity_field_real(c_loc(err), c_null_ptr, c_loc(name_wleb), values)
         call check_api_error(error, err, "Cavity handle is missing")
      end if
      if (.not. allocated(error)) then
         if (any(values /= -12345.0_c_double)) then
            call test_failed(error, "get_cavity_field_real wrote for a missing handle")
         end if
      end if

      deallocate (err)

   end subroutine test_cavity_field_null_handle

   !> A cavity handle with a null inner pointer is reported, not dereferenced
   !>
   !> This is the state a host holds between `moist_new_*_cavity` and the first
   !> `moist_update_cavity`: the handle is real, the cavity behind it is not
   subroutine test_cavity_field_null_inner(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav
      character(kind=c_char), allocatable, target :: name_wleb(:)
      integer(c_int) :: nfield
      real(c_double) :: values(1)

      allocate (err)
      allocate (cav)
      name_wleb = c_string("wleb")

      nfield = -1_c_int
      call moist_get_cavity_field_count(c_loc(err), c_loc(cav), nfield)
      call check_api_error(error, err, "Cavity is not initialized")
      if (.not. allocated(error)) then
         if (nfield /= -1_c_int) then
            call test_failed(error, "get_cavity_field_count changed nfield on failure")
         end if
      end if
      if (allocated(err%ptr)) deallocate (err%ptr)

      if (.not. allocated(error)) then
         values = -12345.0_c_double
         call moist_get_cavity_field_real(c_loc(err), c_loc(cav), c_loc(name_wleb), values)
         call check_api_error(error, err, "Cavity is not initialized")
      end if
      if (.not. allocated(error)) then
         if (any(values /= -12345.0_c_double)) then
            call test_failed(error, "get_cavity_field_real wrote for an empty cavity")
         end if
      end if

      deallocate (cav)
      deallocate (err)

   end subroutine test_cavity_field_null_inner

   !> An empty field name is rejected rather than matched against a declaration
   !>
   !> A bare `""` survives the null-pointer guard and reaches the name scan, so
   !> without its own check it would be compared against every declared name --
   !> harmless today, but only because no cavity declares an empty name
   subroutine test_cavity_field_empty_name(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph
      character(kind=c_char), allocatable, target :: name_empty(:)
      real(c_double) :: values(1)

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      name_empty = c_string("")
      values = -12345.0_c_double
      call moist_get_cavity_field_real(verror, vcav, c_loc(name_empty), values)
      call check_api_error(error, err, "Field name is empty")

      if (.not. allocated(error)) then
         if (any(values /= -12345.0_c_double)) then
            call test_failed(error, "get_cavity_field_real wrote for an empty name")
         end if
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_cavity_field_empty_name

   !> A field index outside the reported count is refused, and the descriptor
   !> outputs remain unchanged
   !>
   !> The index is the one field-API argument the caller derives from an earlier
   !> call, so a stale count is the realistic way to get here.  Unchecked it
   !> would index `query%info` at 0 or past `nfield`; the cleared outputs matter
   !> because a caller that ignores the error must not read a plausible-looking
   !> descriptor out of the buffers
   subroutine test_cavity_field_index_range(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph, nfield
      integer(c_int) :: dtype, rank, dims(2), count
      character(kind=c_char), target :: name(64)

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      call moist_get_cavity_field_count(verror, vcav, nfield)
      if (allocated(err%ptr)) then
         call test_failed(error, "get_cavity_field_count failed: "//err%ptr%message)
      else if (nfield <= 0_c_int) then
         call test_failed(error, "A built cavity declares no fields")
      end if

      if (.not. allocated(error)) then
         call probe_field_info(verror, vcav, -1_c_int, name, dtype, rank, dims, count)
         call check_api_error(error, err, "Field index out of range")
      end if
      if (.not. allocated(error)) then
         call check_unchanged_info(error, "negative index", name, dtype, rank, dims, count)
      end if
      if (allocated(err%ptr)) deallocate (err%ptr)

      if (.not. allocated(error)) then
         call probe_field_info(verror, vcav, nfield, name, dtype, rank, dims, count)
         call check_api_error(error, err, "Field index out of range")
      end if
      if (.not. allocated(error)) then
         call check_unchanged_info(error, "one past the end", name, dtype, rank, dims, count)
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_cavity_field_index_range

   !> Seed the descriptor buffers with detectable values, then describe a field
   subroutine probe_field_info(verror, vcav, ifield, name, dtype, rank, dims, count)
      !> Opaque error and cavity handles
      type(c_ptr), intent(in) :: verror, vcav
      !> Field position to describe
      integer(c_int), intent(in) :: ifield
      !> Name buffer, seeded so any write is visible
      character(kind=c_char), intent(inout) :: name(:)
      !> Descriptor outputs, seeded so a cleared value is distinguishable
      integer(c_int), intent(out) :: dtype, rank, dims(2), count

      name = "Z"
      dtype = -1_c_int
      rank = -1_c_int
      dims = -1_c_int
      count = -1_c_int
      call moist_get_cavity_field_info(verror, vcav, ifield, name, dtype, rank, dims, count)

   end subroutine probe_field_info

   !> Assert a rejected describe left the descriptor and name unchanged
   subroutine check_unchanged_info(error, what, name, dtype, rank, dims, count)
      type(error_type), allocatable, intent(out) :: error
      !> Case being checked, used in the failure message
      character(len=*), intent(in) :: what
      !> Name buffer that must still carry its seed
      character(kind=c_char), intent(in) :: name(:)
      !> Descriptor outputs that must retain their sentinel values
      integer(c_int), intent(in) :: dtype, rank, dims(2), count

      if (dtype /= -1_c_int .or. rank /= -1_c_int .or. count /= -1_c_int) then
         call test_failed(error, "get_cavity_field_info changed the descriptor for "//what)
         return
      end if
      if (any(dims /= -1_c_int)) then
         call test_failed(error, "get_cavity_field_info changed extents for "//what)
         return
      end if
      if (any(name /= "Z")) then
         call test_failed(error, "get_cavity_field_info wrote a name for "//what)
      end if

   end subroutine check_unchanged_info

   !> Description queries report full length while copying only what fits
   subroutine test_cavity_field_about_buffer(error)
      !> Test diagnostic
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph
      character(kind=c_char), allocatable, target :: name_ngrid(:)
      character(kind=c_char), allocatable :: about(:)
      integer(c_size_t) :: nchar, length

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return
      name_ngrid = c_string("ngrid")
      description_checks: block
         call moist_get_cavity_field_about(verror, vcav, c_loc(name_ngrid), &
            & capacity=0_c_size_t, length=nchar)
         call check(error, .not. allocated(err%ptr), "Length-only query failed")
         if (allocated(error)) exit description_checks
         call check(error, nchar > 0, "Expected a nonempty description")
         if (allocated(error)) exit description_checks
         allocate(about(nchar + 1))
         about = "Z"
         call moist_get_cavity_field_about(verror, vcav, c_loc(name_ngrid), about, nchar, length)
         call check(error, .not. allocated(err%ptr), "Truncation should succeed")
         if (allocated(error)) exit description_checks
         call check(error, length == nchar .and. int(c_string_len(about), c_size_t) == nchar - 1)
         if (allocated(error)) exit description_checks
         call check(error, about(nchar + 1) == "Z", "Wrote past supplied capacity")
         if (allocated(error)) exit description_checks
         call moist_get_cavity_field_about(verror, vcav, c_loc(name_ngrid), about, nchar + 1, length)
         call check(error, .not. allocated(err%ptr), "Exact allocation should succeed")
         if (allocated(error)) exit description_checks
         call check(error, length == nchar .and. int(c_string_len(about), c_size_t) == nchar)
      end block description_checks
      call drop_water_cavity(err, vmol, vcav)
   end subroutine test_cavity_field_about_buffer

   !> The logical accessor honours the same contract as the other two
   !>
   !> `converged` is the only boolean a cavity declares, so the bool accessor is
   !> the one branch of `check_field_payload` a C host reaches without a Python
   !> layer in between.  Its payload is cross-checked against the independent
   !> `get_cavity_results` path, since every value it can return is also a legal
   !> seed and a buffer left untouched would otherwise pass unnoticed
   subroutine test_cavity_field_bool(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph, out_ngrid, out_nsph
      real(c_double) :: area, volume
      real(c_double), allocatable :: xyz(:, :), a(:), radii(:), asph(:)
      integer(c_int), allocatable :: owner(:)
      logical(c_bool), allocatable :: reference(:), values(:)
      character(kind=c_char), allocatable, target :: name_conv(:), name_wleb(:)

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      name_conv = c_string("converged")
      name_wleb = c_string("wleb")
      allocate (xyz(3, ngrid), a(ngrid), owner(ngrid), reference(ngrid))
      allocate (radii(nsph), asph(nsph), values(ngrid))

      call moist_get_cavity_results(verror, vcav, ngrid, nsph, area, volume, &
                                    out_ngrid, out_nsph, xyz, a, owner, reference, &
                                    radii, asph)
      if (allocated(err%ptr)) then
         call test_failed(error, "get_cavity_results failed: "//err%ptr%message)
      end if

      !> The declared field has to agree with the fixed getter, point for point
      if (.not. allocated(error)) then
         values = .not. reference
         call moist_get_cavity_field_bool(verror, vcav, c_loc(name_conv), values)
         if (allocated(err%ptr)) then
            call test_failed(error, "get_cavity_field_bool failed: "//err%ptr%message)
         else if (.not. all(values .eqv. reference)) then
            call test_failed(error, "get_cavity_field_bool disagrees with get_cavity_results")
         end if
      end if

      !> Reading a real field through the logical accessor is refused too
      if (.not. allocated(error)) then
         values = .false._c_bool
         call moist_get_cavity_field_bool(verror, vcav, c_loc(name_wleb), values)
         call check_api_error(error, err, "has a different element type")
      end if
      if (.not. allocated(error)) then
         if (any(values)) then
            call test_failed(error, "get_cavity_field_bool wrote for a real-valued field")
         end if
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_cavity_field_bool

   !> Assert that an entry point rejected an undersized capacity, then clear the
   !> error so the handle can be reused
   subroutine expect_capacity_error(error, err, routine)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer, intent(inout) :: err
      character(len=*), intent(in) :: routine

      if (.not. allocated(err%ptr)) then
         call test_failed(error, routine//" accepted an undersized array capacity")
         return
      end if
      if (index(err%ptr%message, "Array capacity") <= 0) then
         call test_failed(error, "Unexpected API error message: "//err%ptr%message)
         return
      end if
      deallocate (err%ptr)

   end subroutine expect_capacity_error

   !> A capacity below the cavity's own ngrid/nsph is refused before a single
   !> element is written.  Without the check these calls write ngrid values into
   !> buffers holding ngrid-1, which is silent heap corruption; the sentinels
   !> below detect any write at all
   subroutine test_capacity_too_small(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph, cap
      integer(c_int) :: out_ngrid, out_nsph
      real(c_double) :: area, volume
      real(c_double), allocatable :: xyz(:, :), a(:), radii(:), asph(:)
      real(c_double), allocatable :: amat0(:, :), xi(:)
      integer(c_int), allocatable :: owner(:)
      logical(c_bool), allocatable :: converged(:)
      real(c_double), parameter :: sentinel = -12345.0_c_double

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      cap = ngrid - 1_c_int
      allocate (xyz(3, cap), a(cap), owner(cap), converged(cap))
      allocate (amat0(cap, cap), xi(cap))
      allocate (radii(nsph), asph(nsph))
      xyz = sentinel
      a = sentinel
      owner = -1_c_int
      converged = .false._c_bool
      amat0 = sentinel
      xi = sentinel
      radii = sentinel
      asph = sentinel

      call moist_get_cavity_results(verror, vcav, cap, nsph, area, volume, &
                                    out_ngrid, out_nsph, xyz, a, owner, converged, &
                                    radii, asph)
      call expect_capacity_error(error, err, "get_cavity_results")
      if (.not. allocated(error)) then
         if (any(a /= sentinel) .or. any(xyz /= sentinel) .or. any(owner /= -1_c_int)) then
            call test_failed(error, "get_cavity_results wrote into a rejected buffer")
         end if
      end if

      if (.not. allocated(error)) then
         call moist_assemble_amat(verror, vcav, cap, amat0, xi)
         call expect_capacity_error(error, err, "assemble_amat")
      end if
      if (.not. allocated(error)) then
         if (any(amat0 /= sentinel) .or. any(xi /= sentinel)) then
            call test_failed(error, "assemble_amat wrote into a rejected buffer")
         end if
      end if

      ! A short per-sphere capacity has to be caught the same way
      if (.not. allocated(error)) then
         deallocate (radii, asph)
         allocate (radii(nsph - 1), asph(nsph - 1))
         radii = sentinel
         asph = sentinel
         call moist_get_cavity_results(verror, vcav, ngrid, nsph - 1_c_int, area, &
                                       volume, out_ngrid, out_nsph, xyz, a, owner, &
                                       converged, radii, asph)
         call expect_capacity_error(error, err, "get_cavity_results (nsph)")
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_capacity_too_small

   !> Accept oversized capacities and preserve padding outside the logical extent
   subroutine test_capacity_oversized(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph, ngrid_cap, nsph_cap
      integer(c_int) :: out_ngrid, out_nsph
      real(c_double) :: area, volume
      real(c_double), allocatable :: xyz(:, :), a(:), radii(:), asph(:)
      integer(c_int), allocatable :: owner(:)
      logical(c_bool), allocatable :: converged(:)
      real(c_double), parameter :: sentinel = -12345.0_c_double

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      ngrid_cap = ngrid + 17_c_int
      nsph_cap = nsph + 5_c_int
      allocate (xyz(3, ngrid_cap), a(ngrid_cap), owner(ngrid_cap), converged(ngrid_cap))
      allocate (radii(nsph_cap), asph(nsph_cap))
      xyz = sentinel
      a = sentinel
      radii = sentinel
      asph = sentinel

      call moist_get_cavity_results(verror, vcav, ngrid_cap, nsph_cap, area, volume, &
                                    out_ngrid, out_nsph, xyz, a, owner, converged, &
                                    radii, asph)

      call check(error, .not. allocated(err%ptr), "Oversized capacity must be accepted")
      if (.not. allocated(error)) call check(error, out_ngrid, ngrid)
      if (.not. allocated(error)) call check(error, out_nsph, nsph)
      if (.not. allocated(error)) call check(error, all(a(:ngrid) /= sentinel))
      if (.not. allocated(error)) call check(error, all(xyz(:, :ngrid) /= sentinel))
      if (.not. allocated(error)) call check(error, all(radii(:nsph) /= sentinel))
      if (.not. allocated(error)) call check(error, all(a(ngrid + 1:) == sentinel))
      if (.not. allocated(error)) call check(error, all(xyz(:, ngrid + 1:) == sentinel))
      if (.not. allocated(error)) call check(error, all(radii(nsph + 1:) == sentinel))

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_capacity_oversized

   !* ================================================================================= *!
   !*                        Isodensity callback failure channel                        *!
   !* ================================================================================= *!

   !> Isodensity LSF callback with a switchable failure channel
   !>
   !> Answers the first `ctx%fail_after` evaluations with a smooth three-Gaussian
   !> model density (so the projection behaves like a real build), then reports
   !> `ctx%fail_status` for every evaluation after that.  The counter is bumped
   !> atomically because moist calls this from inside OpenMP loops
   !>
   !> @param[in]  context Callback context, an `iso_cb_ctx`
   !> @param[in]  point   Evaluation point in Bohr
   !> @param[out] rho     Electron density
   !> @param[out] drho    Density spatial gradient
   !> @param[out] d2rho   Density spatial Hessian (Fortran (3,3), or NULL)
   !> @param[out] d3rho   Density third spatial derivative (Fortran (3,3,3), or NULL)
   !> @returns            0 while healthy, `ctx%fail_status` once failing
   function iso_switchable_callback(context, point, rho, drho, d2rho, d3rho) &
      result(status) bind(C)
      type(c_ptr), value :: context
      real(c_double), intent(in) :: point(3)
      real(c_double), intent(out) :: rho
      real(c_double), intent(out) :: drho(3)
      type(c_ptr), value :: d2rho
      type(c_ptr), value :: d3rho
      integer(c_int) :: status

      type(iso_cb_ctx), pointer :: ctx
      real(c_double), pointer :: hptr(:, :), tptr(:, :, :)
      real(c_double) :: d2rho_l(3, 3), d3rho_l(3, 3, 3)
      real(c_double) :: d(3), g
      integer :: iatom, i, j, k, ncalls
      logical :: want_hess, want_third

      if (.not. c_associated(context)) then
         status = -1_c_int
         return
      end if
      call c_f_pointer(context, ctx)

      !$omp atomic capture
      ctx%calls = ctx%calls + 1
      ncalls = ctx%calls
      !$omp end atomic

      if (ncalls > ctx%fail_after) then
         status = ctx%fail_status
         return
      end if

      status = 0_c_int
      want_hess = c_associated(d2rho)
      want_third = c_associated(d3rho)

      rho = 0.0_c_double
      drho = 0.0_c_double
      d2rho_l = 0.0_c_double
      d3rho_l = 0.0_c_double

      do iatom = 1, 3
         d = point - cb_centers(:, iatom)
         g = cb_c*exp(-cb_a*dot_product(d, d))
         rho = rho + g
         drho = drho - 2.0_c_double*cb_a*d*g
         if (want_hess) then
            do j = 1, 3
               do i = 1, 3
                  d2rho_l(i, j) = d2rho_l(i, j) + (4.0_c_double*cb_a*cb_a*d(i)*d(j) &
                                                   - 2.0_c_double*cb_a*delta(i, j))*g
               end do
            end do
         end if
         if (want_third) then
            do k = 1, 3
               do j = 1, 3
                  do i = 1, 3
                     d3rho_l(i, j, k) = d3rho_l(i, j, k) &
                                        + (-8.0_c_double*cb_a**3*d(i)*d(j)*d(k) &
                                           + 4.0_c_double*cb_a*cb_a*(d(i)*delta(j, k) &
                                                                     + d(j)*delta(i, k) &
                                                                     + d(k)*delta(i, j)))*g
                  end do
               end do
            end do
         end if
      end do

      ! The bare density: moist applies the isovalue and the DROP sign convention
      if (want_hess) then
         call c_f_pointer(d2rho, hptr, [3, 3])
         hptr = d2rho_l
      end if
      if (want_third) then
         call c_f_pointer(d3rho, tptr, [3, 3, 3])
         tptr = d3rho_l
      end if

   end function iso_switchable_callback

   !> Kronecker delta as a c_double, keeping the derivative expressions readable
   pure function delta(i, j) result(d)
      integer, intent(in) :: i, j
      real(c_double) :: d

      d = 0.0_c_double
      if (i == j) d = 1.0_c_double
   end function delta

   !> Build a water cavity backed by `iso_switchable_callback` and update it
   !>
   !> Returns the handles so the caller can inspect the error and delete them
   !>
   !> @param[out]   err        API error handle (caller deallocates)
   !> @param[out]   verror     Opaque view of `err`
   !> @param[out]   vmol,vcav  Structure and cavity handles (caller deletes)
   !> @param[inout] ctx        Callback state; must outlive the cavity handle
   subroutine run_iso_callback_update(err, verror, vmol, vcav, ctx)
      type(vp_error), pointer, intent(out) :: err
      type(c_ptr), intent(out) :: verror
      type(c_ptr), intent(out) :: vmol, vcav
      type(iso_cb_ctx), intent(inout), target :: ctx

      integer(c_int), parameter :: numbers(3) = [8, 1, 1]
      !> Lebedev order, kept small so the build stays cheap
      integer(c_int), target :: nleb

      nleb = 26_c_int
      ctx%calls = 0

      allocate (err)
      verror = c_loc(err)

      vmol = moist_new_structure(verror, 3_c_int, numbers, cb_centers)
      block
         type(c_ptr) :: lsf
         type(iso_options), target :: options
         options = iso_options(c_sizeof(options), cb_rho_iso, 1000.0_c_double)
         lsf = moist_new_isodensity_callback_lsf(verror, c_funloc(iso_switchable_callback), &
                                               c_loc(ctx), c_loc(options))
         vcav = moist_new_drop_cavity(verror, lsf, c_null_ptr, c_null_ptr)
         call moist_delete_lsf(lsf)
      end block
      call moist_update_cavity(verror, vcav, vmol)

   end subroutine run_iso_callback_update

   !> Assert that the update failed with the callback's own status in the message
   !>
   !> @param[out]   error  Test-drive error
   !> @param[in]    err    API error handle produced by the update
   !> @param[in]    needle Status fragment the message has to name
   subroutine expect_iso_callback_error(error, err, needle)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer, intent(in) :: err
      character(len=*), intent(in) :: needle

      if (.not. allocated(err%ptr)) then
         call test_failed(error, "A failing isodensity callback still produced a cavity")
         return
      end if
      if (index(err%ptr%message, "External LSF evaluation failed") <= 0) then
         call test_failed(error, "Unexpected API error message: "//err%ptr%message)
         return
      end if
      if (index(err%ptr%message, needle) <= 0) then
         call test_failed(error, "Error does not name the callback status: "//err%ptr%message)
         return
      end if

   end subroutine expect_iso_callback_error

   !> A callback that fails on its very first evaluation aborts the build
   subroutine test_iso_callback_fails_first(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      type(iso_cb_ctx), target :: ctx

      ctx%fail_after = 0
      ctx%fail_status = 13_c_int

      call run_iso_callback_update(err, verror, vmol, vcav, ctx)
      call expect_iso_callback_error(error, err, "status 13")

      ! The callback is entered once and then not again: the failure latches and
      ! the build unwinds instead of grinding through the rest of the grid
      if (.not. allocated(error)) call check(error, ctx%calls, 1)

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_iso_callback_fails_first

   !> A callback that fails partway through the projection aborts the build
   !>
   !> This is the case that exercises the abort machinery: the projection runs
   !> its anchors in an OpenMP worksharing construct, which cannot be branched
   !> out of, so the failure has to travel out on a shared flag and become an
   !> error only after the region closes.  A first-evaluation failure would pass
   !> even if that were broken
   !>
   !> The callback must also stop being called once it has reported failure, and
   !> a later build with a healthy callback must succeed -- the failure latch is
   !> per-build state, not a permanent poisoning of the cavity
   subroutine test_iso_callback_fails_mid_loop(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      type(iso_cb_ctx), target :: ctx
      integer :: calls_at_abort
      integer(c_int) :: ngrid, nsph
      integer, parameter :: fail_after = 50

      ctx%fail_after = fail_after
      ctx%fail_status = 7_c_int

      call run_iso_callback_update(err, verror, vmol, vcav, ctx)
      calls_at_abort = ctx%calls
      call expect_iso_callback_error(error, err, "status 7")

      if (.not. allocated(error)) then
         if (calls_at_abort <= fail_after) then
            call test_failed(error, "The callback's failing branch was never reached")
         end if
      end if

      ! Once the failure is latched the host callback is not entered again, so
      ! the count stays near the point of failure instead of covering the grid
      if (.not. allocated(error)) then
         if (calls_at_abort > fail_after + 4096) then
            call test_failed(error, "The failing callback kept being called after it reported failure")
         end if
      end if

      ! Recovery: the same cavity handle rebuilds cleanly once the callback
      ! stops failing
      if (.not. allocated(error)) then
         deallocate (err%ptr)
         ctx%calls = 0
         ctx%fail_after = huge(1)
         call moist_update_cavity(verror, vcav, vmol)
         if (allocated(err%ptr)) then
            call test_failed(error, "Cavity stayed broken after the callback recovered: " &
                             //err%ptr%message)
         else
            call moist_get_cavity_sizes(verror, vcav, ngrid, nsph)
            if (.not. allocated(error)) call check(error, nsph, 3_c_int)
            if (.not. allocated(error)) call check(error, ngrid > 0)
         end if
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_iso_callback_fails_mid_loop

   !> Response copies write exactly the logical shape at all ranks, leave the
   !> rest of a larger buffer untouched, and name a missing current item, an
   !> array the current item does not have and a NULL buffer
   !>
   !> @param[out] error Test failure
   subroutine test_response_arrays(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Native wrappers
      type(vp_error), target :: err
      type(vp_response), target :: response
      !> Distinct values expose layout and truncation errors
      type(density_response_type) :: density
      !> Two logical points in four-point buffers
      real(c_double), target :: rho(4), grad(3, 4), hess(3, 3, 4)
      !> NUL-terminated array names
      character(kind=c_char), allocatable, target :: w_rho(:), w_grad_rho(:), &
         & w_hess_rho(:), w_phi(:), bogus(:)
      !> Padding marker
      real(c_double), parameter :: sentinel = -12345.0_c_double
      !> Value index
      integer :: i

      allocate (density%w_rho(2), density%w_grad_rho(3, 2), density%w_hess_rho(3, 3, 2))
      density%w_rho = [1.0_c_double, 2.0_c_double]
      density%w_grad_rho = reshape([(real(i, c_double), i=3, 8)], [3, 2])
      density%w_hess_rho = reshape([(real(i, c_double), i=9, 26)], [3, 3, 2])
      call response_accumulate(response%ptr, density, err%ptr)
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      w_rho = c_string("w_rho")
      w_grad_rho = c_string("w_grad_rho")
      w_hess_rho = c_string("w_hess_rho")
      w_phi = c_string("w_phi")
      bogus = c_string("w_bogus")

      ! Arrays are read from the item the walk stopped at, and none is current yet
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_rho), c_loc(rho))
      call check_api_error(error, err, "[moist_get_response_array] No current response item - call next() first")
      if (allocated(error)) return
      call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
      if (allocated(error)) return

      ! The logical points are written and the rest of the buffer is left alone
      rho = sentinel
      grad = sentinel
      hess = sentinel
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_rho), c_loc(rho))
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_grad_rho), c_loc(grad))
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_hess_rho), c_loc(hess))
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call check(error, all(rho(:2) == density%w_rho))
      if (allocated(error)) return
      call check(error, all(grad(:, :2) == density%w_grad_rho))
      if (allocated(error)) return
      call check(error, all(hess(:, :, :2) == density%w_hess_rho))
      if (allocated(error)) return
      call check(error, all(rho(3:) == sentinel) .and. all(grad(:, 3:) == sentinel) &
         & .and. all(hess(:, :, 3:) == sentinel))
      if (allocated(error)) return

      ! An array the current item does not have and a NULL buffer are named
      call response_get_api(c_loc(err), c_loc(response), c_loc(bogus), c_loc(rho))
      call check_api_error(error, err, "[moist_get_response_array] density has no array 'w_bogus'")
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_phi), c_loc(rho))
      call check_api_error(error, err, "[moist_get_response_array] density has no array 'w_phi'")
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_grad_rho), c_null_ptr)
      call check_api_error(error, err, "Null array pointer provided for 'w_grad_rho'")
      if (allocated(error)) return

      ! The pass ends after the only item, and nothing is current any more
      call check(error, .not. logical(next_response_item_api(c_loc(err), c_loc(response))))
      if (allocated(error)) return
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_rho), c_loc(rho))
      call check_api_error(error, err, "No current response item - call next() first")
   end subroutine test_response_arrays

   !> The C walk visits every item once per pass in accumulation order, rewinds
   !> after the last one, names the current item and fails by name outside a pass
   !>
   !> @param[out] error Test failure
   subroutine test_response_walk(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Native wrappers
      type(vp_error), target :: err
      type(vp_response), target :: response
      !> Items accumulated in this order
      type(potential_adjoint_response_type) :: adjoint
      type(density_response_type) :: density
      !> Pass index
      integer :: pass

      ! A missing response handle ends the loop with an error
      call check(error, .not. logical(next_response_item_api(c_loc(err), c_null_ptr)))
      if (allocated(error)) return
      call check_api_error(error, err, "[moist_next_response_item] Response handle is missing")
      if (allocated(error)) return
      ! A response nothing has filled has nothing to visit, which is no error
      call check(error, .not. logical(next_response_item_api(c_loc(err), c_loc(response))))
      if (allocated(error)) return
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call check(error, current_item_name(err, response), "")
      if (allocated(error)) return
      call check_api_error(error, err, &
         & "[moist_get_response_item_name] No current response item - call next() first")
      if (allocated(error)) return

      adjoint%w_phi = [1.0_c_double, 2.0_c_double]
      allocate (density%w_rho(2), source=3.0_c_double)
      allocate (density%w_grad_rho(3, 2), source=4.0_c_double)
      allocate (density%w_hess_rho(3, 3, 2), source=5.0_c_double)
      call response_accumulate(response%ptr, adjoint, err%ptr)
      call response_accumulate(response%ptr, density, err%ptr)
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return

      ! Two passes: the first ends by rewinding, so the second starts over
      do pass = 1, 2
         call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
         if (allocated(error)) return
         call check(error, current_item_name(err, response), "potential_adjoint")
         if (allocated(error)) return
         call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
         if (allocated(error)) return
         call check(error, current_item_name(err, response), "density")
         if (allocated(error)) return
         call check(error, .not. logical(next_response_item_api(c_loc(err), c_loc(response))), &
            & more="the pass did not end after the last item")
         if (allocated(error)) return
         call check(error, .not. allocated(err%ptr))
         if (allocated(error)) return
      end do

      ! A missing name buffer is refused
      call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
      if (allocated(error)) return
      call response_item_name_api(c_loc(err), c_loc(response))
      call check_api_error(error, err, "[moist_get_response_item_name] Output pointer is missing")
   end subroutine test_response_walk

   !> Name of the current response item through the C entry point, empty on error
   function current_item_name(err, response) result(name)
      !> Error wrapper receiving a failure
      type(vp_error), target, intent(inout) :: err
      !> Response wrapper under test
      type(vp_response), target, intent(inout) :: response
      !> Decoded name
      character(len=:), allocatable :: name
      !> MOIST_NAME_MAX + 1 fits every name
      character(kind=c_char) :: buffer(33)
      integer :: ichar

      buffer = c_null_char
      call response_item_name_api(c_loc(err), c_loc(response), buffer)
      name = ""
      do ichar = 1, c_string_len(buffer)
         name = name//buffer(ichar)
      end do
   end function current_item_name


   !* ================================================================================= *!
   !*                        Validation branches of the C entry points                  *!
   !* ================================================================================= *!

   !> Assert an API error naming `needle`, then clear it so the handle is reused
   subroutine expect_error(error, err, needle)
      !> Test-drive error
      type(error_type), allocatable, intent(out) :: error
      !> API error handle produced by the entry point
      type(vp_error), pointer, intent(in) :: err
      !> Expected message fragment
      character(len=*), intent(in) :: needle

      call check_api_error(error, err, needle)
      if (allocated(err%ptr)) deallocate (err%ptr)

   end subroutine expect_error

   !> Fail when a rejected call changed an output it promises to leave alone
   subroutine check_untouched(error, wrote, what)
      !> Test-drive error, only set when `wrote`
      type(error_type), allocatable, intent(inout) :: error
      !> Whether an output changed
      logical, intent(in) :: wrote
      !> Case being checked, used in the failure message
      character(len=*), intent(in) :: what

      if (allocated(error) .or. .not. wrote) return
      call test_failed(error, "Output changed although the call failed: "//what)

   end subroutine check_untouched

   !> Call one option initializer by position
   subroutine init_options(verror, ikind, buffer)
      !> Opaque error handle, NULL for the guard under test
      type(c_ptr), intent(in) :: verror
      !> Initializer index: drop, iswig, svdw, cfc, isodensity, model, pcm
      integer, intent(in) :: ikind
      !> Caller-owned options storage
      integer(c_int8_t), target, intent(inout) :: buffer(256)

      select case (ikind)
      case (1)
         call moist_init_drop_options(verror, c_loc(buffer), c_sizeof(buffer))
      case (2)
         call moist_init_iswig_options(verror, c_loc(buffer), c_sizeof(buffer))
      case (3)
         call moist_init_svdw_options(verror, c_loc(buffer), c_sizeof(buffer))
      case (4)
         call moist_init_cfc_options(verror, c_loc(buffer), c_sizeof(buffer))
      case (5)
         call moist_init_isodensity_options(verror, c_loc(buffer), c_sizeof(buffer))
      case (6)
         call moist_init_model_options(verror, c_loc(buffer), c_sizeof(buffer))
      case (7)
         call moist_init_pcm_options(verror, c_loc(buffer), c_sizeof(buffer))
      case default
         error stop "test_api: unhandled ikind"
      end select

   end subroutine init_options

   !> Option initializers without an error handle return before touching the buffer
   !>
   !> There is no handle to report into, so the caller's buffer is the only
   !> observable: seeded bytes survive the NULL-handle call, and the same call
   !> with a handle overwrites them, so it is the guard that kept them and not
   !> an initializer that writes nothing
   subroutine test_init_options_null_error(error)
      type(error_type), allocatable, intent(out) :: error
      !> API error handle for the positive control
      type(vp_error), pointer :: err
      !> Caller-owned options storage, larger than every layout
      integer(c_int8_t) :: buffer(256)
      !> Initializer names, in the order of `init_options`
      character(len=10), parameter :: kinds(7) = [character(len=10) :: "drop", "iswig", &
         & "svdw", "cfc", "isodensity", "model", "pcm"]
      integer :: ikind

      allocate (err)
      do ikind = 1, size(kinds)
         buffer = 77_c_int8_t
         call init_options(c_null_ptr, ikind, buffer)
         if (any(buffer /= 77_c_int8_t)) then
            call test_failed(error, "init_"//trim(kinds(ikind))//"_options wrote without an error handle")
            exit
         end if
         call init_options(c_loc(err), ikind, buffer)
         if (allocated(err%ptr)) then
            call test_failed(error, "init_"//trim(kinds(ikind))//"_options failed: "//err%ptr%message)
            exit
         end if
         if (all(buffer == 77_c_int8_t)) then
            call test_failed(error, "init_"//trim(kinds(ikind))//"_options wrote nothing with a handle")
            exit
         end if
      end do
      deallocate (err)

   end subroutine test_init_options_null_error

   !> LSF constructors refuse a missing callback, non-positive isodensity
   !> controls and every malformed basis, and return no handle
   !>
   !> The basis cases walk the shell guard (no shells, a NULL array), both
   !> halves of the primitive guard (a zero count, a sum past `huge`) and the
   !> basis constructor's own rejection of an unsupported angular momentum
   subroutine test_lsf_constructor_guards(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      !> Isodensity options with a negative isovalue
      type(iso_options), target :: options
      !> Two-shell basis, edited per case
      integer(c_int), target :: shell_atom(2), shell_l(2), shell_nprim(2)
      real(c_double), target :: exps(2), coeffs(2)
      !> Returned LSF handle
      type(c_ptr) :: lsf
      !> Expected diagnostics, one per case
      character(len=48), parameter :: needles(7) = [character(len=48) :: &
         & "Callback is missing", "rho_iso and scale must be positive", &
         & "Missing basis arrays or invalid shell count", &
         & "Missing basis arrays or invalid shell count", "Invalid primitive counts", &
         & "Invalid primitive counts", "angular momentum out of supported range"]
      integer :: icase

      allocate (err)
      shell_atom = 0_c_int
      shell_l = 0_c_int
      exps = 0.5_c_double
      coeffs = 1.0_c_double
      do icase = 1, size(needles)
         shell_nprim = 1_c_int
         select case (icase)
         case (1)
            lsf = moist_new_isodensity_callback_lsf(c_loc(err), c_null_funptr, c_null_ptr, c_null_ptr)
         case (2)
            options = iso_options(c_sizeof(options), -1.0_c_double, 1000.0_c_double)
            lsf = moist_new_isodensity_callback_lsf(c_loc(err), c_funloc(iso_switchable_callback), &
                                                  c_null_ptr, c_loc(options))
         case (3)
            lsf = moist_new_isodensity_lsf(c_loc(err), 0_c_int, c_loc(shell_atom), c_loc(shell_l), &
                                           c_loc(shell_nprim), c_loc(exps), c_loc(coeffs), c_null_ptr)
         case (4)
            lsf = moist_new_isodensity_lsf(c_loc(err), 1_c_int, c_loc(shell_atom), c_loc(shell_l), &
                                           c_loc(shell_nprim), c_loc(exps), c_null_ptr, c_null_ptr)
         case (5)
            shell_nprim = [0_c_int, 1_c_int]
            lsf = moist_new_isodensity_lsf(c_loc(err), 2_c_int, c_loc(shell_atom), c_loc(shell_l), &
                                           c_loc(shell_nprim), c_loc(exps), c_loc(coeffs), c_null_ptr)
         case (6)
            ! Each count is valid on its own; only their sum overflows
            shell_nprim = huge(1_c_int)
            lsf = moist_new_isodensity_lsf(c_loc(err), 2_c_int, c_loc(shell_atom), c_loc(shell_l), &
                                           c_loc(shell_nprim), c_loc(exps), c_loc(coeffs), c_null_ptr)
         case (7)
            shell_l = 99_c_int
            lsf = moist_new_isodensity_lsf(c_loc(err), 1_c_int, c_loc(shell_atom), c_loc(shell_l), &
                                           c_loc(shell_nprim), c_loc(exps), c_loc(coeffs), c_null_ptr)
         case default
            error stop "test_api: unhandled icase"
         end select
         call expect_error(error, err, trim(needles(icase)))
         call check_untouched(error, c_associated(lsf), trim(needles(icase)))
         call moist_delete_lsf(lsf)
         if (allocated(error)) exit
      end do
      deallocate (err)

   end subroutine test_lsf_constructor_guards

   !> Cavity constructors refuse each invalid DROP control and a radii handle
   !> whose model was never set, and return no handle
   subroutine test_cavity_constructor_guards(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      !> DROP options, reset to the defaults before each case
      type(drop_options), target :: options
      !> Radii handle that was never given a model
      type(vp_radii), target :: radii
      !> Opaque handles
      type(c_ptr) :: verror, lsf, cavity
      !> Expected diagnostics: six option cases, then the two radii cases
      character(len=48), parameter :: needles(8) = [character(len=48) :: &
         & "Invalid DROP options", "Invalid DROP options", "Invalid DROP options", &
         & "Invalid DROP options", "Invalid DROP options", "Invalid DROP options", &
         & "[moist_new_drop_cavity] Radii are not initialized", &
         & "[moist_new_iswig_cavity] Radii are not initialized"]
      integer :: icase

      allocate (err)
      verror = c_loc(err)
      lsf = moist_new_svdw_lsf(verror, c_null_ptr)
      if (allocated(err%ptr)) then
         call test_failed(error, "LSF setup failed: "//err%ptr%message)
         deallocate (err)
         return
      end if

      do icase = 1, size(needles)
         call moist_init_drop_options(verror, c_loc(options), c_sizeof(options))
         select case (icase)
         case (1)
            options%tolerance = 0.0_c_double
         case (2)
            options%rho_grid_h = 0.0_c_double
         case (3)
            options%branch_weight_s = -1.0_c_double
         case (4)
            options%proj_maxiter = 0_c_int
         case (5)
            options%proj_level = 0_c_int
         case (6)
            options%proj_level = 10_c_int
         case default
            ! Cases 7 and 8 keep the default options
            continue
         end select
         select case (icase)
         case (1:6)
            cavity = moist_new_drop_cavity(verror, lsf, c_null_ptr, c_loc(options))
         case (7)
            cavity = moist_new_drop_cavity(verror, lsf, c_loc(radii), c_null_ptr)
         case (8)
            cavity = moist_new_iswig_cavity(verror, c_loc(radii), c_null_ptr)
         case default
            error stop "test_api: unhandled icase"
         end select
         call expect_error(error, err, trim(needles(icase)))
         call check_untouched(error, c_associated(cavity), trim(needles(icase)))
         call moist_delete_cavity(cavity)
         if (allocated(error)) exit
      end do

      call moist_delete_lsf(lsf)
      deallocate (err)

   end subroutine test_cavity_constructor_guards

   !> Structure and custom-radii entry points name a missing array, a missing
   !> structure and an empty one, and a failed constructor returns no handle
   subroutine test_structure_radii_guards(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      !> Structure that was never constructed
      type(vp_structure), target :: empty
      !> H2, coordinates in Bohr
      integer(c_int) :: numbers(2)
      real(c_double) :: positions(3, 2), radii(2)
      !> Returned structure handle
      type(c_ptr) :: verror, vmol
      !> Expected diagnostics, one per case
      character(len=80), parameter :: needles(8) = [character(len=80) :: &
         & "[moist_new_structure] Required pointer 'numbers' is missing", &
         & "[moist_new_structure] Required pointer 'positions' is missing", &
         & "[moist_update_structure] Required pointer 'positions' is missing", &
         & "[moist_update_structure] Molecular structure data is missing", &
         & "[moist_update_structure] Invalid molecular structure data provided", &
         & "[moist_set_custom_radii_atoms] Required pointer 'atom_radii' is missing", &
         & "[moist_set_custom_radii_elements] Required pointer 'atomic_numbers' is missing", &
         & "[moist_set_custom_radii_elements] Required pointer 'element_radii' is missing"]
      integer :: icase

      numbers = 1_c_int
      positions = reshape([0.0_c_double, 0.0_c_double, 0.0_c_double, &
                           0.0_c_double, 0.0_c_double, 1.4_c_double], [3, 2])
      radii = 2.0_c_double
      allocate (err)
      verror = c_loc(err)
      do icase = 1, size(needles)
         vmol = c_null_ptr
         select case (icase)
         case (1)
            vmol = moist_new_structure(verror, 2_c_int, positions=positions)
         case (2)
            vmol = moist_new_structure(verror, 2_c_int, numbers)
         case (3)
            call update_structure_api(verror, c_null_ptr)
         case (4)
            call update_structure_api(verror, c_null_ptr, positions)
         case (5)
            call update_structure_api(verror, c_loc(empty), positions)
         case (6)
            call set_custom_radii_atoms_api(verror, c_null_ptr, 2_c_int)
         case (7)
            call set_custom_radii_elements_api(verror, c_null_ptr, 2_c_int, element_radii=radii)
         case (8)
            call set_custom_radii_elements_api(verror, c_null_ptr, 2_c_int, numbers)
         case default
            error stop "test_api: unhandled icase"
         end select
         call expect_error(error, err, trim(needles(icase)))
         call check_untouched(error, c_associated(vmol), "structure handle")
         call moist_delete_structure(vmol)
         if (allocated(error)) exit
      end do
      deallocate (err)

   end subroutine test_structure_radii_guards

   !> Omit one required pointer at a time and expect the entry point to name it
   !>
   !> The presence checks precede the handle checks, so the cavity handle is
   !> NULL throughout: a message naming the omitted argument proves its own
   !> branch fired and not the handle guard behind it
   subroutine check_missing_pointers(error, family, names)
      !> Test-drive error
      type(error_type), allocatable, intent(out) :: error
      !> Entry point without the common moist_ prefix
      character(len=*), intent(in) :: family
      !> Required pointer names, in argument order
      character(len=*), intent(in) :: names(:)
      type(vp_error), pointer :: err
      integer :: omit

      allocate (err)
      do omit = 1, size(names)
         call probe_missing_pointer(c_loc(err), family, omit)
         call expect_error(error, err, &
            & "[moist_"//family//"] Required pointer '"//trim(names(omit))//"' is missing")
         if (allocated(error)) exit
      end do
      deallocate (err)

   end subroutine check_missing_pointers

   !> Call `family` with a NULL cavity and its `omit`-th pointer argument absent
   !>
   !> An unallocated allocatable actual argument is absent to an optional
   !> dummy, the Fortran spelling of passing NULL from C; every other buffer is
   !> allocated so only the omitted one can trip a presence check
   subroutine probe_missing_pointer(verror, family, omit)
      !> Opaque error handle
      type(c_ptr), intent(in) :: verror
      !> Entry point without the common moist_ prefix
      character(len=*), intent(in) :: family
      !> Position of the omitted pointer argument
      integer, intent(in) :: omit
      !> Scalar outputs
      real(c_double), allocatable :: area, volume
      integer(c_int), allocatable :: ngrid, nsph, dtype, rank, count
      !> Grid and sphere outputs
      real(c_double), allocatable :: xyz(:, :), a(:), radii(:), asph(:), xi(:), amat0(:, :)
      integer(c_int), allocatable :: owner(:), dims(:)
      logical(c_bool), allocatable :: converged(:)
      character(kind=c_char), allocatable :: name(:)
      !> Gradient outputs
      real(c_double), allocatable :: A_tot1_rA(:, :), V_tot1_rA(:, :), vec1_rA(:, :, :), &
         & xi1_rA(:, :, :), a_i1_rA(:, :, :), v_i1_rA(:, :, :), asph1_rA(:, :, :), &
         & vsph1_rA(:, :, :), r_iI1_rA(:, :, :), xyz1_rA(:, :, :, :), Amat1_rA(:, :, :, :)

      select case (family)
      case ("get_cavity_results")
         allocate (area, volume, ngrid, nsph, xyz(3, 1), a(1), owner(1), converged(1), &
                   radii(1), asph(1))
         if (omit == 1) deallocate (area)
         if (omit == 2) deallocate (volume)
         if (omit == 3) deallocate (ngrid)
         if (omit == 4) deallocate (nsph)
         if (omit == 5) deallocate (xyz)
         if (omit == 6) deallocate (a)
         if (omit == 7) deallocate (owner)
         if (omit == 8) deallocate (converged)
         if (omit == 9) deallocate (radii)
         if (omit == 10) deallocate (asph)
         call moist_get_cavity_results(verror, c_null_ptr, 1_c_int, 1_c_int, area, volume, &
                                       ngrid, nsph, xyz, a, owner, converged, radii, asph)
      case ("get_cavity_sizes")
         allocate (ngrid, nsph)
         if (omit == 1) deallocate (ngrid)
         if (omit == 2) deallocate (nsph)
         call moist_get_cavity_sizes(verror, c_null_ptr, ngrid, nsph)
      case ("get_cavity_field_count")
         call moist_get_cavity_field_count(verror, c_null_ptr)
      case ("get_cavity_field_info")
         allocate (name(65), dtype, rank, dims(2), count)
         if (omit == 1) deallocate (name)
         if (omit == 2) deallocate (dtype)
         if (omit == 3) deallocate (rank)
         if (omit == 4) deallocate (dims)
         if (omit == 5) deallocate (count)
         call moist_get_cavity_field_info(verror, c_null_ptr, 0_c_int, name, dtype, rank, dims, count)
      case ("get_cavity_field_real")
         call moist_get_cavity_field_real(verror, c_null_ptr, c_null_ptr)
      case ("get_cavity_field_int")
         call moist_get_cavity_field_int(verror, c_null_ptr, c_null_ptr)
      case ("get_cavity_field_bool")
         call moist_get_cavity_field_bool(verror, c_null_ptr, c_null_ptr)
      case ("assemble_amat")
         allocate (amat0(1, 1), xi(1))
         if (omit == 1) deallocate (amat0)
         if (omit == 2) deallocate (xi)
         call moist_assemble_amat(verror, c_null_ptr, 1_c_int, amat0, xi)
      case ("get_anchor_gradient")
         allocate (xyz1_rA(3, 3, 1, 1), xi1_rA(3, 1, 1), a_i1_rA(3, 1, 1), v_i1_rA(3, 1, 1), &
                   A_tot1_rA(3, 1), V_tot1_rA(3, 1))
         if (omit == 1) deallocate (xyz1_rA)
         if (omit == 2) deallocate (xi1_rA)
         if (omit == 3) deallocate (a_i1_rA)
         if (omit == 4) deallocate (v_i1_rA)
         if (omit == 5) deallocate (A_tot1_rA)
         if (omit == 6) deallocate (V_tot1_rA)
         call get_anchor_gradient_api(verror, c_null_ptr, 1_c_int, 1_c_int, xyz1_rA, xi1_rA, &
                                      a_i1_rA, v_i1_rA, A_tot1_rA, V_tot1_rA)
      case ("get_cavity_gradient")
         allocate (A_tot1_rA(3, 1), V_tot1_rA(3, 1), asph1_rA(3, 1, 1), vsph1_rA(3, 1, 1), &
                   xyz1_rA(3, 3, 1, 1), r_iI1_rA(3, 1, 1), vec1_rA(3, 1, 1))
         if (omit == 1) deallocate (A_tot1_rA)
         if (omit == 2) deallocate (V_tot1_rA)
         if (omit == 3) deallocate (asph1_rA)
         if (omit == 4) deallocate (vsph1_rA)
         if (omit == 5) deallocate (xyz1_rA)
         if (omit == 6) deallocate (r_iI1_rA)
         if (omit == 7) deallocate (vec1_rA)
         call get_cavity_gradient_api(verror, c_null_ptr, 1_c_int, 1_c_int, A_tot1_rA, V_tot1_rA, &
                                      asph1_rA, vsph1_rA, xyz1_rA, r_iI1_rA, vec1_rA)
      case ("get_amat_gradient")
         allocate (amat0(1, 1), Amat1_rA(3, 1, 1, 1), xi(1))
         if (omit == 1) deallocate (amat0)
         if (omit == 2) deallocate (Amat1_rA)
         if (omit == 3) deallocate (xi)
         call get_amat_gradient_api(verror, c_null_ptr, 1_c_int, 1_c_int, amat0, Amat1_rA, xi)
      case default
         error stop "test_api: unhandled family"
      end select

   end subroutine probe_missing_pointer

   !> Every required output of the result, size, field and A-matrix readers is
   !> named when absent
   subroutine test_missing_pointers_results(error)
      type(error_type), allocatable, intent(out) :: error
      !> Required outputs of each reader, in argument order
      character(len=9), parameter :: results(10) = [character(len=9) :: "area", "volume", &
         & "ngrid", "nsph", "xyz", "a", "owner", "converged", "radii", "asph"]
      character(len=5), parameter :: sizes(2) = [character(len=5) :: "ngrid", "nsph"]
      character(len=6), parameter :: nfield(1) = ["nfield"]
      character(len=5), parameter :: info(5) = [character(len=5) :: "name", "dtype", "rank", &
         & "dims", "count"]
      character(len=6), parameter :: values(1) = ["values"]
      character(len=5), parameter :: amat(2) = [character(len=5) :: "amat0", "xi"]

      call check_missing_pointers(error, "get_cavity_results", results)
      if (allocated(error)) return
      call check_missing_pointers(error, "get_cavity_sizes", sizes)
      if (allocated(error)) return
      call check_missing_pointers(error, "get_cavity_field_count", nfield)
      if (allocated(error)) return
      call check_missing_pointers(error, "get_cavity_field_info", info)
      if (allocated(error)) return
      call check_missing_pointers(error, "get_cavity_field_real", values)
      if (allocated(error)) return
      call check_missing_pointers(error, "get_cavity_field_int", values)
      if (allocated(error)) return
      call check_missing_pointers(error, "get_cavity_field_bool", values)
      if (allocated(error)) return
      call check_missing_pointers(error, "assemble_amat", amat)

   end subroutine test_missing_pointers_results

   !> Every required output of the three gradient readers is named when absent
   subroutine test_missing_pointers_gradients(error)
      type(error_type), allocatable, intent(out) :: error
      !> Required outputs of each reader, in argument order
      character(len=9), parameter :: anchor(6) = [character(len=9) :: "xyz1_rA", "xi1_rA", &
         & "a_i1_rA", "v_i1_rA", "A_tot1_rA", "V_tot1_rA"]
      character(len=9), parameter :: cavity(7) = [character(len=9) :: "A_tot1_rA", "V_tot1_rA", &
         & "asph1_rA", "vsph1_rA", "xyz1_rA", "r_iI1_rA", "rho1_rA"]
      character(len=8), parameter :: amat(3) = [character(len=8) :: "Amat0", "Amat1_rA", "xi"]

      call check_missing_pointers(error, "get_anchor_gradient", anchor)
      if (allocated(error)) return
      call check_missing_pointers(error, "get_cavity_gradient", cavity)
      if (allocated(error)) return
      call check_missing_pointers(error, "get_amat_gradient", amat)

   end subroutine test_missing_pointers_gradients

   !> Create an SvdW DROP cavity and an iSwiG cavity that were never updated
   !>
   !> Both are real handles whose cavities have no surface yet: the state a
   !> host holds between construction and the first `moist_update_cavity`
   subroutine new_unbuilt_cavities(error, err, verror, vdrop, viswig)
      !> Test-drive error
      type(error_type), allocatable, intent(out) :: error
      !> API error handle (allocated here, released by `drop_unbuilt_cavities`)
      type(vp_error), pointer, intent(out) :: err
      !> Opaque view of `err`
      type(c_ptr), intent(out) :: verror
      !> DROP and iSwiG cavity handles
      type(c_ptr), intent(out) :: vdrop, viswig
      type(c_ptr) :: lsf

      allocate (err)
      verror = c_loc(err)
      lsf = moist_new_svdw_lsf(verror, c_null_ptr)
      vdrop = moist_new_drop_cavity(verror, lsf, c_null_ptr, c_null_ptr)
      call moist_delete_lsf(lsf)
      viswig = moist_new_iswig_cavity(verror, c_null_ptr, c_null_ptr)
      if (allocated(err%ptr)) then
         call test_failed(error, "Cavity setup failed: "//err%ptr%message)
         deallocate (err%ptr)
      end if

   end subroutine new_unbuilt_cavities

   !> Release the handles built by new_unbuilt_cavities
   subroutine drop_unbuilt_cavities(err, vdrop, viswig)
      !> API error handle
      type(vp_error), pointer, intent(inout) :: err
      !> DROP and iSwiG cavity handles
      type(c_ptr), intent(inout) :: vdrop, viswig

      call moist_delete_cavity(vdrop)
      call moist_delete_cavity(viswig)
      deallocate (err)

   end subroutine drop_unbuilt_cavities

   !> Readers refuse a cavity that was never built, updates refuse a missing or
   !> empty structure, and the A-matrix assembly refuses a missing or empty
   !> handle; no caller buffer is written
   subroutine test_unbuilt_cavity_guards(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      !> Handle whose cavity pointer is null
      type(vp_cavity), pointer :: empty_cav
      !> Structure that was never constructed
      type(vp_structure), target :: empty_mol
      type(c_ptr) :: verror, vdrop, viswig
      !> Caller buffers, seeded so any write is visible
      real(c_double) :: area, volume, xyz(3, 1), a(1), radii(1), asph(1), amat0(1, 1), xi(1)
      integer(c_int) :: ngrid, nsph, owner(1)
      logical(c_bool) :: converged(1)
      logical :: wrote
      real(c_double), parameter :: sentinel = -12345.0_c_double
      !> Expected diagnostics, one per case
      character(len=64), parameter :: needles(6) = [character(len=64) :: &
         & "[moist_get_cavity_sizes] Cavity is not built yet", &
         & "[moist_get_cavity_results] Cavity is not built yet", &
         & "[moist_update_cavity] Molecular structure data is missing", &
         & "[moist_update_cavity] Invalid number of atoms", &
         & "[moist_assemble_amat] DROP cavity handle is missing", &
         & "[moist_assemble_amat] DROP cavity is not initialized"]
      integer :: icase

      call new_unbuilt_cavities(error, err, verror, vdrop, viswig)
      if (allocated(error)) then
         call drop_unbuilt_cavities(err, vdrop, viswig)
         return
      end if
      allocate (empty_cav)

      do icase = 1, size(needles)
         area = sentinel
         volume = sentinel
         amat0 = sentinel
         xi = sentinel
         ngrid = -1_c_int
         nsph = -1_c_int
         select case (icase)
         case (1)
            call moist_get_cavity_sizes(verror, vdrop, ngrid, nsph)
         case (2)
            call moist_get_cavity_results(verror, vdrop, 1_c_int, 1_c_int, area, volume, ngrid, &
                                          nsph, xyz, a, owner, converged, radii, asph)
         case (3)
            call moist_update_cavity(verror, vdrop, c_null_ptr)
         case (4)
            call moist_update_cavity(verror, vdrop, c_loc(empty_mol))
         case (5)
            call moist_assemble_amat(verror, c_null_ptr, 1_c_int, amat0, xi)
         case (6)
            call moist_assemble_amat(verror, c_loc(empty_cav), 1_c_int, amat0, xi)
         case default
            error stop "test_api: unhandled icase"
         end select
         wrote = area /= sentinel .or. volume /= sentinel .or. ngrid /= -1_c_int &
            & .or. nsph /= -1_c_int .or. any(amat0 /= sentinel) .or. any(xi /= sentinel)
         call expect_error(error, err, trim(needles(icase)))
         call check_untouched(error, wrote, trim(needles(icase)))
         if (allocated(error)) exit
      end do

      deallocate (empty_cav)
      call drop_unbuilt_cavities(err, vdrop, viswig)

   end subroutine test_unbuilt_cavity_guards

   !> DROP-only entry points name a missing handle or output, an empty cavity
   !> and the wrong cavity or LSF type, and leave the tolerance output alone
   subroutine test_drop_only_entry_points(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      !> Handle whose cavity pointer is null
      type(vp_cavity), pointer :: empty_cav
      type(c_ptr) :: verror, vdrop, viswig
      !> Tolerance output and a one-component density
      real(c_double), target :: tolerance, dcart(1, 1)
      real(c_double), parameter :: sentinel = -12345.0_c_double
      !> Expected diagnostics, one per case
      character(len=96), parameter :: needles(6) = [character(len=96) :: &
         & "[moist_get_drop_cavity_tolerance] Cavity handle is missing", &
         & "[moist_get_drop_cavity_tolerance] Output pointer is missing", &
         & "[moist_get_drop_cavity_tolerance] Cavity is not initialized", &
         & "[moist_get_drop_cavity_tolerance] Supplied cavity is not a DROP cavity", &
         & "[moist_set_isodensity_density] Cavity LSF is not the internal isodensity type", &
         & "[moist_set_isodensity_density] Cavity is not DROP type"]
      integer :: icase

      call new_unbuilt_cavities(error, err, verror, vdrop, viswig)
      if (allocated(error)) then
         call drop_unbuilt_cavities(err, vdrop, viswig)
         return
      end if
      allocate (empty_cav)
      dcart = 1.0_c_double

      do icase = 1, size(needles)
         tolerance = sentinel
         select case (icase)
         case (1)
            call get_drop_cavity_tolerance_api(verror, c_null_ptr, c_loc(tolerance))
         case (2)
            call get_drop_cavity_tolerance_api(verror, vdrop, c_null_ptr)
         case (3)
            call get_drop_cavity_tolerance_api(verror, c_loc(empty_cav), c_loc(tolerance))
         case (4)
            call get_drop_cavity_tolerance_api(verror, viswig, c_loc(tolerance))
         case (5)
            call set_isodensity_density_api(verror, vdrop, 1_c_int, c_loc(dcart))
         case (6)
            call set_isodensity_density_api(verror, viswig, 1_c_int, c_loc(dcart))
         case default
            error stop "test_api: unhandled icase"
         end select
         call expect_error(error, err, trim(needles(icase)))
         call check_untouched(error, tolerance /= sentinel, trim(needles(icase)))
         if (allocated(error)) exit
      end do

      ! The same reader answers for a DROP cavity, so the refusals above are
      ! the type guard and not a reader that never writes
      if (.not. allocated(error)) then
         call get_drop_cavity_tolerance_api(verror, vdrop, c_loc(tolerance))
         if (allocated(err%ptr)) then
            call test_failed(error, "get_drop_cavity_tolerance failed: "//err%ptr%message)
         else if (.not. tolerance > 0.0_c_double) then
            call test_failed(error, "get_drop_cavity_tolerance returned no positive tolerance")
         end if
      end if

      deallocate (empty_cav)
      call drop_unbuilt_cavities(err, vdrop, viswig)

   end subroutine test_drop_only_entry_points

   !> The three surface contractions name a missing or empty handle, NULL
   !> buffers, and a cavity without the surface data they contract
   !>
   !> An unbuilt DROP cavity has zero grid points, so the one-element buffers
   !> below are never indexed; they only have to be non-NULL
   subroutine test_contraction_guards(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      !> Handle whose cavity pointer is null
      type(vp_cavity), pointer :: empty_cav
      type(c_ptr) :: verror, vdrop, viswig
      !> Non-NULL stand-ins for every array argument
      real(c_double), target :: w1(1), w2(1), w3(3), w4(1), w5(3), w6(9)
      !> Expected diagnostics, one per case
      character(len=104), parameter :: needles(9) = [character(len=104) :: &
         & "[moist_contract_surface_lsf_weights_extended] Cavity is not initialized", &
         & "[moist_contract_surface_lsf_weights_extended] Null array pointer provided", &
         & "[moist_contract_surface_lsf_weights_extended] contract_surface_lsf_weights: "// &
         & "cavity surface data are incomplete", &
         & "[moist_contract_pcm_nuclear_gradient] Cavity handle is missing", &
         & "[moist_contract_pcm_nuclear_gradient] Cavity is not initialized", &
         & "[moist_contract_pcm_nuclear_gradient] Cavity position derivatives are unavailable", &
         & "[moist_contract_amat1_q1q2_surface_weights] Cavity handle is missing", &
         & "[moist_contract_amat1_q1q2_surface_weights] Cavity is not initialized", &
         & "[moist_contract_amat1_q1q2_surface_weights] Cavity does not provide a Gaussian PCM surface"]
      integer :: icase

      call new_unbuilt_cavities(error, err, verror, vdrop, viswig)
      if (allocated(error)) then
         call drop_unbuilt_cavities(err, vdrop, viswig)
         return
      end if
      allocate (empty_cav)

      do icase = 1, size(needles)
         select case (icase)
         case (1)
            call contract_surface_lsf_weights_extended_api(verror, c_loc(empty_cav), c_null_ptr, &
               & c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, &
               & c_null_ptr, c_null_ptr)
         case (2)
            call contract_surface_lsf_weights_extended_api(verror, vdrop, c_null_ptr, &
               & c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, &
               & c_null_ptr, c_null_ptr)
         case (3)
            call contract_surface_lsf_weights_extended_api(verror, vdrop, c_loc(w1), c_loc(w2), &
               & c_loc(w3), c_loc(w4), c_loc(w5), c_loc(w6), c_null_ptr, c_null_ptr, c_null_ptr)
         case (4)
            call contract_pcm_nuclear_gradient_api(verror, c_null_ptr, c_null_ptr, c_null_ptr, &
               & c_null_ptr, c_null_ptr)
         case (5)
            call contract_pcm_nuclear_gradient_api(verror, c_loc(empty_cav), c_null_ptr, &
               & c_null_ptr, c_null_ptr, c_null_ptr)
         case (6)
            call contract_pcm_nuclear_gradient_api(verror, vdrop, c_loc(w1), c_loc(w3), &
               & c_loc(w2), c_loc(w5))
         case (7)
            call contract_amat1_q1q2_surface_weights_api(verror, c_null_ptr, c_null_ptr, &
               & c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr)
         case (8)
            call contract_amat1_q1q2_surface_weights_api(verror, c_loc(empty_cav), c_null_ptr, &
               & c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr)
         case (9)
            call contract_amat1_q1q2_surface_weights_api(verror, vdrop, c_loc(w1), c_loc(w2), &
               & c_loc(w4), c_loc(w1), c_loc(w3))
         case default
            error stop "test_api: unhandled icase"
         end select
         call expect_error(error, err, trim(needles(icase)))
         if (allocated(error)) exit
      end do

      deallocate (empty_cav)
      call drop_unbuilt_cavities(err, vdrop, viswig)

   end subroutine test_contraction_guards

   !> Wrap a hand-filled stub cavity in a cavity handle
   !>
   !> Two grid points on one sphere; the first width is zero, the position
   !> derivatives are sized for two spheres, and the only declared field has a
   !> name one character longer than MOIST_FIELD_NAME_MAX
   subroutine new_stub_cavity(cav)
      !> Handle owning the stub (release with deallocate, not moist_delete_cavity)
      type(vp_cavity), pointer, intent(out) :: cav

      allocate (cav)
      allocate (stub_cavity :: cav%ptr)
      select type (stub => cav%ptr)
      type is (stub_cavity)
         stub%field_name = repeat("x", 65)
         stub%ngrid = 2
         stub%nsph = 1
         stub%xi0 = [0.0_wp, 1.0_wp]
         stub%f = [1.0_wp, 1.0_wp]
         stub%xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp], [3, 2])
         allocate (stub%sphxyz(3, 1), source=0.0_wp)
         allocate (stub%xyz1_rA(3, 3, 2, 2), source=0.0_wp)
      end select

   end subroutine new_stub_cavity

   !> Refusals only a hand-built cavity reaches: an overlong field name, a
   !> non-positive width in the A-matrix kernels, and position derivatives
   !> whose extents disagree with the sphere count
   subroutine test_stub_cavity_guards(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav
      type(c_ptr) :: verror, vstub
      !> Field descriptor outputs
      integer(c_int) :: dtype, rank, dims(2), count
      character(kind=c_char) :: name(66)
      !> Caller buffers, seeded so any write is visible
      real(c_double), target :: amat0(2, 2), xi(2), q(2), w_xi(2), w_f(2), w_xyz(3, 2)
      real(c_double), target :: za(1), grad(3, 1)
      logical :: wrote
      real(c_double), parameter :: sentinel = -12345.0_c_double
      !> Expected diagnostics, one per case
      character(len=104), parameter :: needles(6) = [character(len=104) :: &
         & "[moist_get_cavity_field_info] Field name exceeds MOIST_FIELD_NAME_MAX", &
         & "[moist_assemble_amat] Gaussian PCM widths must be positive", &
         & "[moist_contract_amat1_q1q2_surface_weights] Null array pointer provided", &
         & "[moist_contract_amat1_q1q2_surface_weights] Gaussian PCM widths must be positive", &
         & "[moist_contract_pcm_nuclear_gradient] Null array pointer provided", &
         & "[moist_contract_pcm_nuclear_gradient] pcm_electrostatic_nuclear_gradient: "// &
         & "array shape mismatch"]
      integer :: icase

      allocate (err)
      verror = c_loc(err)
      call new_stub_cavity(cav)
      vstub = c_loc(cav)
      q = 1.0_c_double
      za = 1.0_c_double

      do icase = 1, size(needles)
         amat0 = sentinel
         xi = sentinel
         w_f = sentinel
         grad = sentinel
         wrote = .false.
         select case (icase)
         case (1)
            call probe_field_info(verror, vstub, 0_c_int, name, dtype, rank, dims, count)
         case (2)
            call moist_assemble_amat(verror, vstub, 2_c_int, amat0, xi)
            wrote = any(amat0 /= sentinel) .or. any(xi /= sentinel)
         case (3)
            call contract_amat1_q1q2_surface_weights_api(verror, vstub, c_loc(q), c_loc(q), &
               & c_null_ptr, c_loc(w_f), c_loc(w_xyz))
            wrote = any(w_f /= sentinel)
         case (4)
            call contract_amat1_q1q2_surface_weights_api(verror, vstub, c_loc(q), c_loc(q), &
               & c_loc(w_xi), c_loc(w_f), c_loc(w_xyz))
         case (5)
            call contract_pcm_nuclear_gradient_api(verror, vstub, c_loc(q), c_null_ptr, &
               & c_loc(za), c_loc(grad))
            wrote = any(grad /= sentinel)
         case (6)
            call contract_pcm_nuclear_gradient_api(verror, vstub, c_loc(q), c_loc(w_xyz), &
               & c_loc(za), c_loc(grad))
         case default
            error stop "test_api: unhandled icase"
         end select
         call expect_error(error, err, trim(needles(icase)))
         call check_untouched(error, wrote, trim(needles(icase)))
         if (icase == 1 .and. .not. allocated(error)) then
            call check_unchanged_info(error, "an overlong name", name, dtype, rank, dims, count)
         end if
         if (allocated(error)) exit
      end do

      deallocate (cav%ptr)
      deallocate (cav)
      deallocate (err)

   end subroutine test_stub_cavity_guards

   !> Model entry points name a missing or empty cavity, a missing or empty
   !> model, a model that is not the general one, a model not updated yet and
   !> an invalid density buffer, and a failed constructor returns no handle
   subroutine test_model_handle_guards(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      !> Handle whose cavity pointer is null
      type(vp_cavity), pointer :: empty_cav
      !> Model handles that were never constructed and one holding a stub
      type(vp_model), pointer :: empty_model, stub
      type(c_ptr) :: verror, vdrop, viswig, vmodel, handle
      !> One-component density
      real(c_double), target :: density(1)
      !> Expected diagnostics, one per case
      character(len=80), parameter :: needles(12) = [character(len=80) :: &
         & "[moist_new_model] Cavity handle is missing", &
         & "[moist_new_model] Cavity is not initialized", &
         & "[moist_update_model] Model is not initialized", &
         & "[moist_update_model] Molecular structure data is missing", &
         & "[moist_get_model_cavity] Model handle is missing", &
         & "[moist_get_model_cavity] Model is not initialized", &
         & "[moist_get_model_cavity] This solvation model type does not expose a cavity", &
         & "[moist_new_coupling] Model is not initialized", &
         & "[moist_new_coupling] Model is not a general solvation model", &
         & "[moist_new_coupling] General model must be updated first", &
         & "[moist_set_model_isodensity_density] Invalid density buffer", &
         & "[moist_set_model_isodensity_density] Invalid density buffer"]
      integer :: icase

      call new_unbuilt_cavities(error, err, verror, vdrop, viswig)
      if (allocated(error)) then
         call drop_unbuilt_cavities(err, vdrop, viswig)
         return
      end if
      vmodel = moist_new_model(verror, viswig, c_null_ptr)
      if (allocated(err%ptr)) then
         call test_failed(error, "Model setup failed: "//err%ptr%message)
         call drop_unbuilt_cavities(err, vdrop, viswig)
         return
      end if
      allocate (empty_cav, empty_model, stub)
      allocate (stub_model :: stub%ptr)
      density = 1.0_c_double

      do icase = 1, size(needles)
         handle = c_null_ptr
         select case (icase)
         case (1)
            handle = moist_new_model(verror, c_null_ptr, c_null_ptr)
         case (2)
            handle = moist_new_model(verror, c_loc(empty_cav), c_null_ptr)
         case (3)
            call update_solvation_model_api(verror, c_loc(empty_model), c_null_ptr)
         case (4)
            call update_solvation_model_api(verror, vmodel, c_null_ptr)
         case (5)
            handle = get_solvation_model_cavity_api(verror, c_null_ptr)
         case (6)
            handle = get_solvation_model_cavity_api(verror, c_loc(empty_model))
         case (7)
            handle = get_solvation_model_cavity_api(verror, c_loc(stub))
         case (8)
            handle = new_coupling_api(verror, c_loc(empty_model))
         case (9)
            handle = new_coupling_api(verror, c_loc(stub))
         case (10)
            handle = new_coupling_api(verror, vmodel)
         case (11)
            call moist_set_model_isodensity_density(verror, vmodel, 0_c_int, c_loc(density))
         case (12)
            call moist_set_model_isodensity_density(verror, vmodel, 1_c_int, c_null_ptr)
         case default
            error stop "test_api: unhandled icase"
         end select
         call expect_error(error, err, trim(needles(icase)))
         call check_untouched(error, c_associated(handle), trim(needles(icase)))
         if (allocated(error)) exit
      end do

      deallocate (stub%ptr)
      deallocate (stub, empty_model, empty_cav)
      call delete_solvation_model_api(vmodel)
      call drop_unbuilt_cavities(err, vdrop, viswig)

   end subroutine test_model_handle_guards

   !> Refusals that need an updated general model: a component added after the
   !> update, a gradient read on a model that was never updated, and a NULL
   !> gradient buffer
   !>
   !> Also pins the NULL-error guard of `moist_set_model_isodensity_density`:
   !> with a handle the call would invalidate the model, so a coupling that
   !> still builds afterwards shows it returned before touching anything
   subroutine test_model_phase_guards(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, viswig, vmodel, vfresh, vpv, vcpl, vresp
      !> H2, coordinates in Bohr
      integer(c_int) :: numbers(2)
      real(c_double) :: positions(3, 2)
      !> Gradient buffer, seeded so any write is visible, and a density
      real(c_double), target :: gradient(3, 2), density(1)
      real(c_double), parameter :: sentinel = -12345.0_c_double

      numbers = 1_c_int
      positions = reshape([0.0_c_double, 0.0_c_double, 0.0_c_double, &
                           0.0_c_double, 0.0_c_double, 1.4_c_double], [3, 2])
      density = 1.0_c_double
      allocate (err)
      verror = c_loc(err)
      vcpl = c_null_ptr
      vresp = c_null_ptr
      vfresh = c_null_ptr

      vmol = moist_new_structure(verror, 2_c_int, numbers, positions)
      viswig = moist_new_iswig_cavity(verror, c_null_ptr, c_null_ptr)
      vmodel = moist_new_model(verror, viswig, c_null_ptr)
      vpv = new_pv_component_api(verror, 1.0e-4_c_double)
      call general_model_add_component_api(verror, vmodel, vpv)
      call update_solvation_model_api(verror, vmodel, vmol)

      checks: block
         if (allocated(err%ptr)) then
            call test_failed(error, "Model setup failed: "//err%ptr%message)
            exit checks
         end if

         call general_model_add_component_api(verror, vmodel, vpv)
         call expect_error(error, err, &
            & "[moist_add_model_component] Components cannot be added after model update")
         if (allocated(error)) exit checks

         call moist_set_model_isodensity_density(c_null_ptr, vmodel, 1_c_int, c_loc(density))
         vcpl = new_coupling_api(verror, vmodel)
         if (allocated(err%ptr)) then
            call test_failed(error, "Coupling on the updated model failed: "//err%ptr%message)
            exit checks
         end if

         vfresh = moist_new_model(verror, viswig, c_null_ptr)
         vresp = new_response_api(verror)
         gradient = sentinel
         call general_model_get_gradient_api(verror, vfresh, vcpl, vresp, 2_c_int, c_loc(gradient))
         call expect_error(error, err, "[moist_get_model_gradient] General model must be updated first")
         call check_untouched(error, any(gradient /= sentinel), "gradient of a model never updated")
         if (allocated(error)) exit checks

         call general_model_get_gradient_api(verror, vmodel, vcpl, vresp, 2_c_int, c_null_ptr)
         call expect_error(error, err, "[moist_get_model_gradient] Null gradient pointer provided")
      end block checks

      call delete_coupling_api(vcpl)
      call delete_response_api(vresp)
      call delete_solvation_model_api(vfresh)
      call delete_solvation_model_api(vmodel)
      call delete_solvation_component_api(vpv)
      call moist_delete_cavity(viswig)
      call moist_delete_structure(vmol)
      deallocate (err)

   end subroutine test_model_phase_guards

   !> Response reads name an empty or unterminated array name and an array the
   !> current item was accumulated without, at every rank, and leave the
   !> caller's buffer alone
   subroutine test_response_unavailable_arrays(error)
      type(error_type), allocatable, intent(out) :: error
      !> Native wrappers
      type(vp_error), target :: err
      type(vp_response), target :: response
      !> Items with only part of their arrays
      type(potential_adjoint_response_type) :: adjoint
      type(density_response_type) :: density
      !> Output buffer, seeded so any write is visible
      real(c_double), target :: values(3, 3, 2)
      !> NUL-terminated array names
      character(kind=c_char), allocatable, target :: empty(:), overlong(:), w_phi(:), &
         & w_grad_rho(:), w_hess_rho(:)
      real(c_double), parameter :: sentinel = -12345.0_c_double

      allocate (density%w_rho(2), source=1.0_c_double)
      call response_accumulate(response%ptr, adjoint, err%ptr)
      call response_accumulate(response%ptr, density, err%ptr)
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      empty = c_string("")
      overlong = c_string(repeat("w", 40))
      w_phi = c_string("w_phi")
      w_grad_rho = c_string("w_grad_rho")
      w_hess_rho = c_string("w_hess_rho")
      values = sentinel

      checks: block
         ! The name is decoded before the walk is consulted, so no item is needed
         call response_get_api(c_loc(err), c_loc(response), c_loc(empty), c_loc(values))
         call expect_error(error, err, "[moist_get_response_array] Invalid output name")
         if (allocated(error)) exit checks
         call response_get_api(c_loc(err), c_loc(response), c_loc(overlong), c_loc(values))
         call expect_error(error, err, "[moist_get_response_array] Invalid output name")
         if (allocated(error)) exit checks

         call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
         if (allocated(error)) exit checks
         call response_get_api(c_loc(err), c_loc(response), c_loc(w_phi), c_loc(values))
         call expect_error(error, err, "[moist_get_response_array] 'w_phi' is not available")
         if (allocated(error)) exit checks

         call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
         if (allocated(error)) exit checks
         call response_get_api(c_loc(err), c_loc(response), c_loc(w_grad_rho), c_loc(values))
         call expect_error(error, err, "[moist_get_response_array] 'w_grad_rho' is not available")
         if (allocated(error)) exit checks
         call response_get_api(c_loc(err), c_loc(response), c_loc(w_hess_rho), c_loc(values))
         call expect_error(error, err, "[moist_get_response_array] 'w_hess_rho' is not available")
         if (allocated(error)) exit checks

         call check_untouched(error, any(values /= sentinel), "response array buffer")
      end block checks

   end subroutine test_response_unavailable_arrays

   !> Gradient readers refuse capacities one short of the built cavity before
   !> writing a single element
   !>
   !> The sphere capacity is short for the anchor reader and the grid capacity
   !> for the other two, so both halves of the capacity checks are exercised;
   !> the PCM nuclear contraction then runs to completion on the same cavity
   subroutine test_gradient_capacity_too_small(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph
      !> Undersized caller buffers
      real(c_double), allocatable :: xyz1_rA(:, :, :, :), xi1_rA(:, :, :), a_i1_rA(:, :, :), &
         & v_i1_rA(:, :, :), A_tot1_rA(:, :), V_tot1_rA(:, :), asph1_rA(:, :, :), &
         & vsph1_rA(:, :, :), r_iI1_rA(:, :, :), rho1_rA(:, :, :), amat0(:, :), &
         & Amat1_rA(:, :, :, :), xi(:)
      !> Zero surface weights, unit charges and the contracted gradient
      real(c_double), allocatable, target :: w_phi(:), w_xyz(:, :), za(:), grad(:, :)
      real(c_double), parameter :: sentinel = -12345.0_c_double

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) then
         call drop_water_cavity(err, vmol, vcav)
         return
      end if

      checks: block
         call compute_anchor_gradient_api(verror, vcav)
         if (allocated(err%ptr)) then
            call test_failed(error, "compute_anchor_gradient failed: "//err%ptr%message)
            exit checks
         end if
         allocate (xyz1_rA(3, 3, nsph - 1, ngrid), xi1_rA(3, nsph - 1, ngrid), &
                   a_i1_rA(3, nsph - 1, ngrid), v_i1_rA(3, nsph - 1, ngrid), &
                   A_tot1_rA(3, nsph - 1), V_tot1_rA(3, nsph - 1), source=sentinel)
         call get_anchor_gradient_api(verror, vcav, nsph - 1_c_int, ngrid, xyz1_rA, xi1_rA, &
                                      a_i1_rA, v_i1_rA, A_tot1_rA, V_tot1_rA)
         call expect_error(error, err, "[moist_get_anchor_gradient] Array capacity is too small")
         call check_untouched(error, any(xyz1_rA /= sentinel) .or. any(xi1_rA /= sentinel) &
            & .or. any(A_tot1_rA /= sentinel), "anchor gradient buffers")
         if (allocated(error)) exit checks
         deallocate (xyz1_rA, A_tot1_rA, V_tot1_rA)

         call compute_cavity_gradient_api(verror, vcav)
         if (allocated(err%ptr)) then
            call test_failed(error, "compute_cavity_gradient failed: "//err%ptr%message)
            exit checks
         end if
         allocate (A_tot1_rA(3, nsph), V_tot1_rA(3, nsph), asph1_rA(3, nsph, nsph), &
                   vsph1_rA(3, nsph, nsph), xyz1_rA(3, 3, nsph, ngrid - 1), &
                   r_iI1_rA(3, nsph, ngrid - 1), rho1_rA(3, nsph, ngrid - 1), source=sentinel)
         call get_cavity_gradient_api(verror, vcav, nsph, ngrid - 1_c_int, A_tot1_rA, V_tot1_rA, &
                                      asph1_rA, vsph1_rA, xyz1_rA, r_iI1_rA, rho1_rA)
         call expect_error(error, err, "[moist_get_cavity_gradient] Array capacity is too small")
         call check_untouched(error, any(A_tot1_rA /= sentinel) .or. any(xyz1_rA /= sentinel) &
            & .or. any(rho1_rA /= sentinel), "cavity gradient buffers")
         if (allocated(error)) exit checks

         allocate (amat0(ngrid - 1, ngrid - 1), Amat1_rA(3, nsph, ngrid - 1, ngrid - 1), &
                   xi(ngrid - 1), source=sentinel)
         call get_amat_gradient_api(verror, vcav, nsph, ngrid - 1_c_int, amat0, Amat1_rA, xi)
         call expect_error(error, err, "[moist_get_amat_gradient] Array capacity is too small")
         call check_untouched(error, any(amat0 /= sentinel) .or. any(Amat1_rA /= sentinel) &
            & .or. any(xi /= sentinel), "A-matrix gradient buffers")
         if (allocated(error)) exit checks

         ! The PCM contraction runs to completion on the same cavity: zero
         ! potential adjoint and position weights leave exactly nothing
         allocate (w_phi(ngrid), w_xyz(3, ngrid), source=0.0_c_double)
         allocate (za(nsph), source=1.0_c_double)
         allocate (grad(3, nsph), source=sentinel)
         call contract_pcm_nuclear_gradient_api(verror, vcav, c_loc(w_phi), c_loc(w_xyz), &
            & c_loc(za), c_loc(grad))
         if (allocated(err%ptr)) then
            call test_failed(error, "contract_pcm_nuclear_gradient failed: "//err%ptr%message)
         else if (any(grad /= 0.0_c_double)) then
            call test_failed(error, "Zero weights gave a nonzero PCM nuclear gradient")
         end if
      end block checks

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_gradient_capacity_too_small

   !> The A-matrix surface weights of a built cavity carry the closed-form
   !> switching-factor channel and no net position force, a NULL output is
   !> refused before the others are written, and a later success clears the
   !> stale error
   !>
   !> - only the diagonal A_ii = sqrt(2/pi)*xi_i/f_i depends on f, so
   !>   w_f_i = -q1_i*q2_i*sqrt(2/pi)*xi_i/f_i**2 exactly
   !> - A depends on point separations only, so a rigid shift of every point
   !>   leaves q1^T A q2 unchanged and the position weights sum to zero
   subroutine test_amat_surface_weights(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph
      !> Surface data read back by name
      real(c_double), allocatable :: xi0(:), f(:)
      !> Contraction vectors and weight outputs
      real(c_double), allocatable, target :: q1(:), q2(:), w_xi(:), w_f(:), w_xyz(:, :)
      !> Closed-form switching-factor weights
      real(c_double), allocatable :: reference(:)
      character(kind=c_char), allocatable, target :: name_xi0(:), name_f(:)
      real(c_double), parameter :: sentinel = -12345.0_c_double
      real(c_double), parameter :: sqrt_2_over_pi = 0.7978845608028654_c_double
      integer :: i

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) then
         call drop_water_cavity(err, vmol, vcav)
         return
      end if
      name_xi0 = c_string("xi0")
      name_f = c_string("f")
      allocate (xi0(ngrid), f(ngrid), q1(ngrid), q2(ngrid), w_xi(ngrid), w_f(ngrid), &
                w_xyz(3, ngrid))

      checks: block
         call moist_get_cavity_field_real(verror, vcav, c_loc(name_xi0), xi0)
         call moist_get_cavity_field_real(verror, vcav, c_loc(name_f), f)
         if (allocated(err%ptr)) then
            call test_failed(error, "Reading the Gaussian surface failed: "//err%ptr%message)
            exit checks
         end if
         ! Distinct vectors, so the q1_i*q2_j + q1_j*q2_i symmetrisation is exercised
         do i = 1, ngrid
            q1(i) = (0.1_c_double + 0.1_c_double*real(i - 1, c_double)/real(ngrid - 1, c_double))*f(i)
            q2(i) = cos(real(i - 1, c_double))*f(i)
         end do

         w_xi = sentinel
         w_f = sentinel
         call contract_amat1_q1q2_surface_weights_api(verror, vcav, c_loc(q1), c_loc(q2), &
            & c_loc(w_xi), c_loc(w_f), c_null_ptr)
         call check_api_error(error, err, &
            & "[moist_contract_amat1_q1q2_surface_weights] Null array pointer provided")
         call check_untouched(error, any(w_xi /= sentinel) .or. any(w_f /= sentinel), &
            & "surface weights with a NULL position output")
         if (allocated(error)) exit checks

         ! The stale error from the refusal above is left in place on purpose
         call contract_amat1_q1q2_surface_weights_api(verror, vcav, c_loc(q1), c_loc(q2), &
            & c_loc(w_xi), c_loc(w_f), c_loc(w_xyz))
         if (allocated(err%ptr)) then
            call test_failed(error, "contract_amat1_q1q2_surface_weights failed: "//err%ptr%message)
            exit checks
         end if

         reference = -q1*q2*sqrt_2_over_pi*xi0/f**2
         if (any(abs(w_f - reference) > 1.0e-12_c_double*abs(reference))) then
            call test_failed(error, "Switching-factor weights disagree with the closed form")
            exit checks
         end if
         if (any(w_xi == sentinel)) then
            call test_failed(error, "Width weights were not written")
            exit checks
         end if
         if (maxval(abs(sum(w_xyz, dim=2))) > 1.0e-10_c_double*maxval(abs(w_xyz))) then
            call test_failed(error, "Position weights carry a net force on a rigid shift")
            exit checks
         end if
      end block checks

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_amat_surface_weights

   !* ================================================================================= *!
   !*                              Test double procedures                               *!
   !* ================================================================================= *!

   !> Stub update: the test sets every array by hand
   subroutine stub_cavity_update(self, mol, error)
      !> Stub cavity
      class(stub_cavity), intent(inout) :: self
      !> Ignored structure
      type(structure_type), intent(in) :: mol
      !> Never set
      type(moist_error_type), allocatable, intent(out) :: error
   end subroutine stub_cavity_update

   !> Stub gradient: the test sets the derivative arrays by hand
   subroutine stub_cavity_gradient(self, error)
      !> Stub cavity
      class(stub_cavity), intent(inout) :: self
      !> Never set
      type(moist_error_type), allocatable, intent(out) :: error
   end subroutine stub_cavity_gradient

   !> Declare the widths under the stub's own field name
   subroutine stub_cavity_fields(self, query)
      !> Stub cavity
      class(stub_cavity), intent(in) :: self
      !> Walker collecting or fetching the declarations
      type(cavity_field_query_type), intent(inout) :: query

      call query%add_real(self%field_name, "Stub widths", self%xi0)

   end subroutine stub_cavity_fields

   !> Stub update: accept any structure
   subroutine stub_model_update(self, mol, error)
      !> Stub model
      class(stub_model), intent(inout) :: self
      !> Ignored structure
      class(structure_type), intent(in) :: mol
      !> Never set
      type(moist_error_type), allocatable, intent(out) :: error
   end subroutine stub_model_update

   !> Stub energy: contribute nothing
   subroutine stub_model_energy(self, coupling, energy, error)
      !> Stub model
      class(stub_model), intent(inout) :: self
      !> Ignored coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Unchanged energy
      real(wp), intent(inout) :: energy
      !> Never set
      type(moist_error_type), allocatable, intent(out) :: error
   end subroutine stub_model_energy

   !> Stub response: contribute nothing
   subroutine stub_model_response(self, coupling, response, error)
      !> Stub model
      class(stub_model), intent(inout) :: self
      !> Ignored coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Unchanged response
      type(response_type), intent(inout) :: response
      !> Never set
      type(moist_error_type), allocatable, intent(out) :: error
   end subroutine stub_model_response

   !> Stub gradient: contribute nothing
   subroutine stub_model_gradient(self, coupling, response, gradient, error)
      !> Stub model
      class(stub_model), intent(inout) :: self
      !> Ignored coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Unchanged response
      type(response_type), intent(inout) :: response
      !> Unchanged gradient
      real(wp), intent(inout) :: gradient(:, :)
      !> Never set
      type(moist_error_type), allocatable, intent(out) :: error
   end subroutine stub_model_gradient

end module test_api
