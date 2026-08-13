! TODO - FUTURE ENHANCEMENTS:
! ---------------------------
! 1. PARAMETER CONFIGURATION:
!    Currently cavity parameters (num_leb, rho_param, switching functions, etc.)
!    are hardcoded or use defaults. Need API to configure parameters:
!    - Option A: Pass params in constructor via opaque parameter handle
!    - Option B: Add moist_set_cavity_params() functions
!    - Option C: Use JSON/TOML configuration file path in constructor
!
! 2. ADDITIONAL CAVITY TYPES:
!    - Implement moist_new_numsa_cavity() constructor
!    - Add corresponding type-specific getters if needed
!
! 3. GRADIENT API:
!    - Expose get_gradient() deferred procedure through C API
!    - Generic moist_get_cavity_gradient() for all types
!
! 4. ADVANCED FEATURES:
!    - Cavity serialization/deserialization
!    - Cavity visualization data export
!    - Performance profiling hooks
!===============================================================================!

module moist_api
   use iso_c_binding
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_structure, only: structure_type, new
   use moist_type, only: cavity_type, solvation_model_type, solvation_model_component_type
   use moist_channels, only: coupling_type, response_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, &
                                             assemble_pcm_amat_with_gradient, pcm_amat_surface_weights, &
                                             pcm_amat_nuclear_gradient
   use moist_model_component_pcm_electrostatics, only: &
      pcm_electrostatic_nuclear_gradient
   use moist_model_component_pcm_type, only: potential_source
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
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_isodensity_callback, only: &
      moist_cavity_drop_lsf_isodensity_callback_type
   use moist_cavity_drop_lsf_isodensity_internal, only: &
      moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_data_solvents, only: solvation_system_parameters, new_solvation_system_parameters, get_solvent_id
   use moist_version, only: get_moist_version
   use moist_output_ascii, only: ascii_moist_header => moist_header, &
      & HEADER_FULL, HEADER_SHORT, HEADER_ASCII, &
      & ascii_moist_build_header => moist_build_header
   implicit none
   private

   character(len=*), parameter :: namespace = "moist_"
   integer, parameter :: api_max_cstr = 4096

   type :: vp_error
      type(error_type), allocatable :: ptr
   end type vp_error

   type :: vp_structure
      type(structure_type) :: ptr
   end type vp_structure

   type :: vp_cavity
      !> Run context owned by this handle; the cavity borrows a pointer to it.
      type(moist_context_type) :: ctx
      class(cavity_type), pointer :: ptr => null()
      logical :: owned = .true.
   end type vp_cavity

   type :: vp_radii
      class(radius_type), allocatable :: ptr
   end type vp_radii

   type :: vp_model
      !> Run context owned by this handle; the model borrows a pointer to it.
      !> A model constructor calls `new_context` on it and hands it to the
      !> concrete model (and on to its components), the way the cavity handles
      !> do below. It is torn down in `delete_solvation_model_api`.
      type(moist_context_type) :: ctx
      class(solvation_model_type), allocatable :: ptr
      !> Host coupling data owned by this handle (traces, response arrays)
      type(coupling_type) :: coupling
   end type vp_model

   type :: vp_component
      !> Run context owned by this handle until the component is copied into a
      !> model; `general_model_add_component` re-points the copy at the model
      !> context so the copy stays valid after this handle is deleted.
      type(moist_context_type) :: ctx
      !> Concrete component owned by this opaque handle
      class(solvation_model_component_type), allocatable :: ptr
   end type vp_component

   public :: vp_error, vp_structure, vp_cavity, vp_radii, vp_model, vp_component
   public :: get_version_api, get_version_string_api
   public :: new_error_api, check_error_api, get_error_api, delete_error_api
   public :: new_structure_api, delete_structure_api, update_structure_api
   public :: new_cpcm_radii_api, new_smd_radii_api, new_d3_radii_api, new_cosmo_radii_api, new_bondi_radii_api
   public :: new_custom_radii_api, set_custom_radii_atoms_api, set_custom_radii_elements_api
   public :: delete_radii_api
   ! Solvation model API
   public :: update_solvation_model_api
   public :: get_solvation_model_energy_api
   public :: get_solvation_model_cavity_api
   public :: delete_solvation_model_api
   ! General solvation model and its components
   public :: new_cpcm_component_api, new_cosmo_component_api
   public :: new_pv_component_api, new_gostshyp_component_api
   public :: delete_solvation_component_api
   public :: new_general_solvation_model_api, general_model_add_component_api
   public :: general_model_supply_electrostatics_api, general_model_supply_gostshyp_api
   public :: general_model_get_trace_response_api, general_model_get_response_api
   public :: general_model_get_response_extended_api
   public :: general_model_get_gradient_api
   ! Type-specific constructors
   public :: new_drop_cavity_api
   public :: new_drop_cavity_with_radii_api
   public :: new_cfc_drop_cavity_api
   public :: new_iswig_cavity_api
   public :: new_drop_cavity_isodensity_callback_api
   public :: new_drop_cavity_isodensity_internal_api
   ! Generic cavity operations
   public :: update_cavity_api
   public :: get_cavity_sizes_api
   public :: get_cavity_results_api
   public :: delete_cavity_api
   ! Type-specific getters
   public :: get_drop_specific_api
   public :: get_drop_numbering_api
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
   public :: contract_nuc_elec_qefield_rA_api
   public :: print_header_api, print_header_api_short, print_header_api_ascii, print_version_api
   public :: print_build_header_api

contains

!> API error helper - creates consistent error messages with routine context
   subroutine api_error(error, routine, msg)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: routine
      character(len=*), intent(in) :: msg
      call fatal_error(error, "["//routine//"] "//msg)
   end subroutine api_error

!> Report a response channel the caller asked for but this model cannot produce
!>
!> A non-null output pointer is a requirement, not a hint: the caller is told
!> that nothing feeds the channel rather than receiving a zero-filled buffer it
!> cannot tell apart from a genuine zero. Passing null skips the channel.
   subroutine api_channel_absent(error, routine, channel)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: routine
      character(len=*), intent(in) :: channel
      call api_error(error, routine, "Requested response channel '"//channel// &
                     & "' is not produced by this model configuration; pass NULL to skip it")
   end subroutine api_channel_absent

!> Obtain library version as major * 10000 + minor * 100 + patch
   function get_version_api() result(version) &
         & bind(C, name=namespace//"get_version")
      integer(c_int) :: version
      integer :: major, minor, patch

      call get_moist_version(major, minor, patch)
      version = 10000_c_int*major + 100_c_int*minor + patch

   end function get_version_api

!> Get version string (e.g., "0.1.0")
   subroutine get_version_string_api(charptr, buffersize) &
         & bind(C, name=namespace//"get_version_string")
      character(kind=c_char), intent(inout) :: charptr(*)
      integer(c_int), intent(in), optional :: buffersize
      integer :: major, minor, patch, max_length
      character(len=32) :: version_str

      if (present(buffersize)) then
         max_length = max(1, buffersize)
      else
         ! Without a C-side bound, only write the terminating null byte.
         max_length = 1
      end if

      call get_moist_version(major, minor, patch)
      write (version_str, '(i0,".",i0,".",i0)') major, minor, patch
      call f_c_character(trim(version_str), charptr, max_length)

   end subroutine get_version_string_api

!> Print MOIST header banner to a file unit (use 6 for stdout)
   subroutine print_header_api(unit) &
         & bind(C, name=namespace//"print_header")
      integer(c_int), value, intent(in) :: unit
      call ascii_moist_header(unit, HEADER_FULL)
   end subroutine print_header_api

!> Print MOIST header short version to a file unit (use 6 for stdout)
   subroutine print_header_api_short(unit) &
         & bind(C, name=namespace//"print_header_short")
      integer(c_int), value, intent(in) :: unit
      call ascii_moist_header(unit, HEADER_SHORT)
   end subroutine print_header_api_short

!> Print MOIST ASCII banner to a file unit (use 6 for stdout)
   subroutine print_header_api_ascii(unit) &
         & bind(C, name=namespace//"print_header_ascii")
      integer(c_int), value, intent(in) :: unit
      call ascii_moist_header(unit, HEADER_ASCII)
   end subroutine print_header_api_ascii

!> Print MOIST version string to a file unit
   subroutine print_version_api(unit) &
         & bind(C, name=namespace//"print_version")
      integer(c_int), value, intent(in) :: unit
      integer :: major, minor, patch

      call get_moist_version(major, minor, patch)
      write (unit, "(a,1x,i0,a,i0,a,i0)") "moist version", major, ".", minor, ".", patch

   end subroutine print_version_api

!> Print moist build information (version, git commit, compiler, host) to a file unit
   subroutine print_build_header_api(unit) &
         & bind(C, name=namespace//"print_build_header")
      integer(c_int), value, intent(in) :: unit
      call ascii_moist_build_header(unit)
   end subroutine print_build_header_api

   !> Create new error handle object
   function new_error_api() &
         & result(verror) &
         & bind(C, name=namespace//"new_error")
      type(vp_error), pointer :: error
      type(c_ptr) :: verror

      allocate (error)
      verror = c_loc(error)

   end function new_error_api

!> Delete error handle object
   subroutine delete_error_api(verror) &
         & bind(C, name=namespace//"delete_error")
      type(c_ptr), intent(inout) :: verror
      type(vp_error), pointer :: error

      if (c_associated(verror)) then
         call c_f_pointer(verror, error)

         deallocate (error)
         verror = c_null_ptr
      end if

   end subroutine delete_error_api

!> Check error handle status
   function check_error_api(verror) result(status) &
         & bind(C, name=namespace//"check_error")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      integer(c_int) :: status

      if (c_associated(verror)) then
         call c_f_pointer(verror, error)

         if (allocated(error%ptr)) then
            status = 1
         else
            status = 0
         end if
      else
         status = 2
      end if

   end function check_error_api

!> Get error message from error handle
   subroutine get_error_api(verror, charptr, buffersize) &
         & bind(C, name=namespace//"get_error")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      character(kind=c_char), intent(inout) :: charptr(*)
      integer(c_int), intent(in), optional :: buffersize
      integer :: max_length

      if (c_associated(verror)) then
         call c_f_pointer(verror, error)

         if (present(buffersize)) then
            max_length = max(1, buffersize)
         else
            ! Without a C-side bound, only write the terminating null byte.
            max_length = 1
         end if

         if (allocated(error%ptr)) then
            call f_c_character(error%ptr%message, charptr, max_length)
         end if
      end if

   end subroutine get_error_api

!> Create new molecular structure data (quantities in Bohr)
   function new_structure_api(verror, natoms, numbers, positions, &
         & c_lattice, c_periodic) result(vmol) &
         & bind(C, name=namespace//"new_structure")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      integer(c_int), value, intent(in) :: natoms
      integer(c_int), intent(in) :: numbers(natoms)
      real(c_double), intent(in) :: positions(3, natoms)
      real(c_double), intent(in), optional :: c_lattice(3, 3)
      real(wp), allocatable :: lattice(:, :)
      logical(c_bool), intent(in), optional :: c_periodic(3)
      logical, allocatable :: periodic(:)
      type(vp_structure), pointer :: mol
      type(c_ptr) :: vmol

      vmol = c_null_ptr

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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

   end function new_structure_api

!> Delete molecular structure data
   subroutine delete_structure_api(vmol) &
         & bind(C, name=namespace//"delete_structure")
      type(c_ptr), intent(inout) :: vmol
      type(vp_structure), pointer :: mol

      if (c_associated(vmol)) then
         call c_f_pointer(vmol, mol)

         deallocate (mol)
         vmol = c_null_ptr
      end if

   end subroutine delete_structure_api

!> Update coordinates and lattice parameters (quantities in Bohr)
   subroutine update_structure_api(verror, vmol, positions, lattice) &
         & bind(C, name=namespace//"update_structure")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vmol
      type(vp_structure), pointer :: mol
      real(c_double), intent(in) :: positions(3, *)
      real(c_double), intent(in), optional :: lattice(3, 3)

      if (.not. c_associated(verror)) then
         return
      end if
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmol)) then
         call api_error(error%ptr, "update_structure_api", "Molecular structure data is missing")
         return
      end if
      call c_f_pointer(vmol, mol)

      if (mol%ptr%nat <= 0 .or. mol%ptr%nid <= 0 .or. .not. allocated(mol%ptr%num) &
         & .or. .not. allocated(mol%ptr%id) .or. .not. allocated(mol%ptr%xyz)) then
         call api_error(error%ptr, "update_structure_api", "Invalid molecular structure data provided")
         return
      end if

      mol%ptr%xyz(:, :) = positions(:3, :mol%ptr%nat)
      if (present(lattice)) then
         mol%ptr%lattice(:, :) = lattice(:3, :3)
      end if

      call verify_structure(error%ptr, mol%ptr)

   end subroutine update_structure_api

!> Create a new radii handle from a named model.
   subroutine new_radii_handle_api(verror, model_name, routine_name, vradii)
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      character(len=*), intent(in) :: model_name
      character(len=*), intent(in) :: routine_name
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

!> Create CPCM radii model handle.
   function new_cpcm_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_cpcm_radii")
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "cpcm", "new_cpcm_radii_api", vradii)
   end function new_cpcm_radii_api

!> Create SMD radii model handle.
   function new_smd_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_smd_radii")
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "smd", "new_smd_radii_api", vradii)
   end function new_smd_radii_api

!> Create D3 radii model handle.
   function new_d3_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_d3_radii")
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "d3", "new_d3_radii_api", vradii)
   end function new_d3_radii_api

!> Create COSMO radii model handle.
   function new_cosmo_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_cosmo_radii")
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "cosmo", "new_cosmo_radii_api", vradii)
   end function new_cosmo_radii_api

!> Create Bondi radii model handle.
   function new_bondi_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_bondi_radii")
      type(c_ptr), value :: verror
      type(c_ptr) :: vradii
      call new_radii_handle_api(verror, "bondi", "new_bondi_radii_api", vradii)
   end function new_bondi_radii_api

!> Create custom radii model handle.
!> Must be initialized with set_custom_radii_atoms or set_custom_radii_elements before use.
   function new_custom_radii_api(verror) result(vradii) &
         & bind(C, name=namespace//"new_custom_radii")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(vp_radii), pointer :: radii
      type(c_ptr) :: vradii

      vradii = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      allocate (radii)
      allocate (radius_type_custom :: radii%ptr)
      vradii = c_loc(radii)
   end function new_custom_radii_api

!> Set custom radii from per-atom values (bohr).
!> @param[in] verror      Error handle
!> @param[in] vradii      Custom radii handle
!> @param[in] natoms      Number of atoms in atom_radii
!> @param[in] atom_radii  Per-atom radii values in bohr
   subroutine set_custom_radii_atoms_api(verror, vradii, natoms, atom_radii) &
         & bind(C, name=namespace//"set_custom_radii_atoms")
      type(c_ptr), value :: verror
      type(c_ptr), value :: vradii
      integer(c_int), value :: natoms
      real(c_double), intent(in) :: atom_radii(*)
      type(vp_error), pointer :: error
      type(vp_radii), pointer :: radii
      type(error_type), allocatable :: radii_error
      real(wp), allocatable :: atom_radii_wp(:)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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

!> Set custom radii from per-element values (bohr).
!> @param[in] verror          Error handle
!> @param[in] vradii          Custom radii handle
!> @param[in] nentries        Number of entries in atomic_numbers and element_radii
!> @param[in] atomic_numbers  Atomic numbers for supplied radii
!> @param[in] element_radii   Per-element radii values in bohr
   subroutine set_custom_radii_elements_api(verror, vradii, nentries, atomic_numbers, element_radii) &
         & bind(C, name=namespace//"set_custom_radii_elements")
      type(c_ptr), value :: verror
      type(c_ptr), value :: vradii
      integer(c_int), value :: nentries
      integer(c_int), intent(in) :: atomic_numbers(*)
      real(c_double), intent(in) :: element_radii(*)
      type(vp_error), pointer :: error
      type(vp_radii), pointer :: radii
      type(error_type), allocatable :: radii_error
      integer, allocatable :: atomic_numbers_f(:)
      real(wp), allocatable :: element_radii_wp(:)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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

!> Delete radii model handle.
   subroutine delete_radii_api(vradii) &
         & bind(C, name=namespace//"delete_radii")
      type(c_ptr), intent(inout) :: vradii
      type(vp_radii), pointer :: radii

      if (c_associated(vradii)) then
         call c_f_pointer(vradii, radii)
         if (allocated(radii%ptr)) deallocate (radii%ptr)
         deallocate (radii)
         vradii = c_null_ptr
      end if
   end subroutine delete_radii_api

!> Delete solvation model handle.
   subroutine delete_solvation_model_api(vmodel) &
         & bind(C, name=namespace//"delete_solvation_model")
      type(c_ptr), intent(inout) :: vmodel
      type(vp_model), pointer :: model

      if (c_associated(vmodel)) then
         call c_f_pointer(vmodel, model)

         if (allocated(model%ptr)) deallocate (model%ptr)
         call model%ctx%delete()
         deallocate (model)
         vmodel = c_null_ptr
      end if

   end subroutine delete_solvation_model_api

!> Update a solvation model with a molecular structure.
   subroutine update_solvation_model_api(verror, vmodel, vmol) &
         & bind(C, name=namespace//"update_solvation_model")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vmodel
      type(vp_model), pointer :: model
      type(c_ptr), value :: vmol
      type(vp_structure), pointer :: mol
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmodel)) then
         call api_error(error%ptr, "update_solvation_model_api", "Model handle is missing")
         return
      end if
      call c_f_pointer(vmodel, model)

      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "update_solvation_model_api", "Model is not initialized")
         return
      end if

      if (.not. c_associated(vmol)) then
         call api_error(error%ptr, "update_solvation_model_api", "Molecular structure data is missing")
         return
      end if
      call c_f_pointer(vmol, mol)

      !> Host coupling data belongs to the geometry and surface it was built on.
      !> A new update must require every external trace and response array to be
      !> supplied again even when the point count happens to stay unchanged.
      call model%coupling%clear()

      call model%ptr%update(mol%ptr, error=model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "update_solvation_model_api", model_error%message)
      end if

   end subroutine update_solvation_model_api

!> Get the total solvation energy from a solvation model
   subroutine get_solvation_model_energy_api(verror, vmodel, energy) &
         & bind(C, name=namespace//"get_solvation_model_energy")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vmodel
      type(vp_model), pointer :: model
      real(c_double), intent(out) :: energy
      type(error_type), allocatable :: model_error

      energy = 0.0_c_double

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmodel)) then
         call api_error(error%ptr, "get_solvation_model_energy_api", "Model handle is missing")
         return
      end if
      call c_f_pointer(vmodel, model)

      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "get_solvation_model_energy_api", "Model is not initialized")
         return
      end if

      call model%ptr%get_energy(model%coupling, energy, model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "get_solvation_model_energy_api", model_error%message)
      end if

   end subroutine get_solvation_model_energy_api

!> Get a borrowed cavity handle from a solvation model.
!> The returned handle is NOT owned by the caller. Independent cavity-update
!> entry points reject it; moist_delete_cavity releases only the borrowed wrapper.
!> It remains valid as long as the parent model exists.
   function get_solvation_model_cavity_api(verror, vmodel) result(vcav) &
         & bind(C, name=namespace//"get_solvation_model_cavity")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vmodel
      type(vp_model), pointer :: model
      type(c_ptr) :: vcav
      type(vp_cavity), pointer :: cav
      class(cavity_type), pointer :: cavity_ptr
      character(len=:), allocatable :: message

      vcav = c_null_ptr

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmodel)) then
         call api_error(error%ptr, "get_solvation_model_cavity", "Model handle is missing")
         return
      end if
      call c_f_pointer(vmodel, model)

      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "get_solvation_model_cavity", "Model is not initialized")
         return
      end if

      call borrow_general_cavity(model%ptr, cavity_ptr, message)
      if (allocated(message)) then
         call api_error(error%ptr, "get_solvation_model_cavity", message)
         return
      end if

      allocate (cav)
      cav%ptr => cavity_ptr
      cav%owned = .false.
      vcav = c_loc(cav)

   end function get_solvation_model_cavity_api

   !> Point at the cavity a general solvation model owns, without taking it
   subroutine borrow_general_cavity(model, cavity_ptr, message)
      !> Solvation model that may own a cavity. `intent(inout)` because the
      !> borrowed pointer stays live and the handle it feeds is read-write.
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

!> Allocate either PCM-family component behind the common opaque handle.
   subroutine new_pcm_component_common(verror, epsilon, solver, use_cosmo, &
                                       routine_name, vcomponent)
      !> Error handle
      type(c_ptr), value :: verror
      !> Relative dielectric constant
      real(c_double), value :: epsilon
      !> PCM solver enumeration
      integer(c_int), value :: solver
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
            call new_component_cosmo(pcm, component%ctx, real(epsilon, wp), &
                                     solver=int(solver), phi_source=potential_source%external, &
                                     error=component_error)
         end select
      else
         allocate (solvation_model_component_cpcm :: item)
         select type (pcm => item)
         type is (solvation_model_component_cpcm)
            call new_component_cpcm(pcm, component%ctx, real(epsilon, wp), &
                                    solver=int(solver), phi_source=potential_source%external, &
                                    error=component_error)
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

!> Create a CPCM component handle for use with a general model.
   function new_cpcm_component_api(verror, epsilon, solver) result(vcomponent) &
         & bind(C, name=namespace//"new_cpcm_component")
      type(c_ptr), value :: verror
      real(c_double), value :: epsilon
      integer(c_int), value :: solver
      type(c_ptr) :: vcomponent

      call new_pcm_component_common(verror, epsilon, solver, .false., &
                                    "new_cpcm_component", vcomponent)

   end function new_cpcm_component_api

!> Create a COSMO component handle for use with a general model.
   function new_cosmo_component_api(verror, epsilon, solver) result(vcomponent) &
         & bind(C, name=namespace//"new_cosmo_component")
      type(c_ptr), value :: verror
      real(c_double), value :: epsilon
      integer(c_int), value :: solver
      type(c_ptr) :: vcomponent

      call new_pcm_component_common(verror, epsilon, solver, .true., &
                                    "new_cosmo_component", vcomponent)

   end function new_cosmo_component_api

!> Create a pressure-volume energy component handle.
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

      allocate (component)
      call new_context(component%ctx, verbosity=0, debug=.false.)
      call new_component_pv(item, real(pressure, wp))
      item%ctx => component%ctx
      allocate (component%ptr, source=item)
      vcomponent = c_loc(component)

   end function new_pv_component_api

!> Create a GOSTSHYP hydrostatic-pressure component handle.
!>
!> The component cannot form its own density traces; the host must supply the
!> Gaussian moments through `general_model_supply_gostshyp` after every cavity
!> update, and read the amplitudes back with
!> `general_model_get_potential_extended`.
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

      allocate (component)
      call new_context(component%ctx, verbosity=0, debug=.false.)
      call new_component_gostshyp(item, real(pressure, wp))
      item%ctx => component%ctx
      allocate (component%ptr, source=item)
      vcomponent = c_loc(component)

   end function new_gostshyp_component_api

!> Delete a solvation-component handle.
   subroutine delete_solvation_component_api(vcomponent) &
         & bind(C, name=namespace//"delete_solvation_component")
      !> Component handle
      type(c_ptr), intent(inout) :: vcomponent
      !> Decoded component wrapper
      type(vp_component), pointer :: component

      if (c_associated(vcomponent)) then
         call c_f_pointer(vcomponent, component)
         if (allocated(component%ptr)) deallocate (component%ptr)
         call component%ctx%delete()
         deallocate (component)
         vcomponent = c_null_ptr
      end if

   end subroutine delete_solvation_component_api

!> Create a general solvation model around an owned copy of a cavity.
   function new_general_solvation_model_api(verror, vcavity, c_debug, c_verbose) result(vmodel) &
         & bind(C, name=namespace//"new_general_solvation_model")
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
         call api_error(error%ptr, "new_general_solvation_model", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcavity, cavity)
      if (.not. associated(cavity%ptr)) then
         call api_error(error%ptr, "new_general_solvation_model", "Cavity is not initialized")
         return
      end if

      allocate (model)
      call new_context(model%ctx, verbosity=int(c_verbose), debug=logical(c_debug))
      call new_model_general(general, cavity%ptr, model%ctx, model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "new_general_solvation_model", model_error%message)
         call model%ctx%delete()
         deallocate (model)
         return
      end if
      allocate (model%ptr, source=general)
      vmodel = c_loc(model)

   end function new_general_solvation_model_api

!> Append a component to a general model.
   subroutine general_model_add_component_api(verror, vmodel, vcomponent) &
         & bind(C, name=namespace//"general_model_add_component")
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

      if (.not. c_associated(vmodel) .or. .not. c_associated(vcomponent)) then
         call api_error(error%ptr, "general_model_add_component", &
                        "Model or component handle is missing")
         return
      end if
      call c_f_pointer(vmodel, model)
      call c_f_pointer(vcomponent, component)
      if (.not. allocated(model%ptr) .or. .not. allocated(component%ptr)) then
         call api_error(error%ptr, "general_model_add_component", &
                        "Model or component is not initialized")
         return
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         ! The model stores a copy, so hand it the model-owned context first:
         ! the copy keeps whatever pointer the source carried, and this handle's
         ! own context dies with `moist_delete_solvation_component`.
         component%ptr%ctx => model%ctx
         call general%add_component(component%ptr, model_error)
         component%ptr%ctx => component%ctx
         if (allocated(model_error)) then
            call api_error(error%ptr, "general_model_add_component", model_error%message)
         end if
      class default
         call api_error(error%ptr, "general_model_add_component", &
                        "Model is not a general solvation model")
      end select

   end subroutine general_model_add_component_api

!> Supply electrostatic traces and response arrays to a general model.
   subroutine general_model_supply_electrostatics_api(verror, vmodel, ngrid, c_phi, &
                                                      c_w_xi, c_w_f, c_w_xyz, c_w_n, c_qefield) &
      & bind(C, name=namespace//"general_model_supply_electrostatics")
      !> Error and model handles
      type(c_ptr), value :: verror, vmodel
      !> Electrostatic grid size
      integer(c_int), value :: ngrid
      !> Input arrays
      type(c_ptr), value :: c_phi, c_w_xi, c_w_f, c_w_xyz, c_w_n, c_qefield
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded model wrapper
      type(vp_model), pointer :: model
      !> Scalar input views
      real(c_double), pointer :: phi(:), w_xi(:), w_f(:)
      !> Vector input views
      real(c_double), pointer :: w_xyz(:, :), w_n(:, :), qefield(:, :)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmodel) .or. .not. c_associated(c_phi) .or. ngrid <= 0) then
         call api_error(error%ptr, "general_model_supply_electrostatics", &
                        "Invalid model, potential, or grid size")
         return
      end if
      call c_f_pointer(vmodel, model)
      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "general_model_supply_electrostatics", &
                        "Model is not initialized")
         return
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         if (.not. general%updated .or. ngrid /= general%cavity%ngrid) then
            call api_error(error%ptr, "general_model_supply_electrostatics", &
                           "Grid size does not match the updated general-model cavity")
            return
         end if
      class default
         call api_error(error%ptr, "general_model_supply_electrostatics", &
                        "Model is not a general solvation model")
         return
      end select

      call c_f_pointer(c_phi, phi, [int(ngrid)])
      model%coupling%electrostatics%phi = real(phi, wp)
      if (allocated(model%coupling%electrostatics%w_xi)) then
        deallocate (model%coupling%electrostatics%w_xi)
      end if
      if (allocated(model%coupling%electrostatics%w_f)) then
        deallocate (model%coupling%electrostatics%w_f)
      end if
      if (allocated(model%coupling%electrostatics%w_xyz)) then
        deallocate (model%coupling%electrostatics%w_xyz)
      end if
      if (allocated(model%coupling%electrostatics%w_normal)) then
        deallocate (model%coupling%electrostatics%w_normal)
      end if
      if (allocated(model%coupling%electrostatics%qefield)) then
        deallocate (model%coupling%electrostatics%qefield)
      end if
      if (c_associated(c_w_xi)) then
         call c_f_pointer(c_w_xi, w_xi, [int(ngrid)])
         model%coupling%electrostatics%w_xi = real(w_xi, wp)
      end if
      if (c_associated(c_w_f)) then
         call c_f_pointer(c_w_f, w_f, [int(ngrid)])
         model%coupling%electrostatics%w_f = real(w_f, wp)
      end if
      if (c_associated(c_w_xyz)) then
         call c_f_pointer(c_w_xyz, w_xyz, [3, int(ngrid)])
         model%coupling%electrostatics%w_xyz = real(w_xyz, wp)
      end if
      if (c_associated(c_w_n)) then
         call c_f_pointer(c_w_n, w_n, [3, int(ngrid)])
         model%coupling%electrostatics%w_normal = real(w_n, wp)
      end if
      if (c_associated(c_qefield)) then
         call c_f_pointer(c_qefield, qefield, [3, int(ngrid)])
         model%coupling%electrostatics%qefield = real(qefield, wp)
      end if

   end subroutine general_model_supply_electrostatics_api

!> Supply the Gaussian density moments the GOSTSHYP component consumes.
!>
!> `gt`, `pt`, `mt` and `rt` are the moments of the solute density against the
!> unnormalized Gaussian `exp(-w_i |r - r_i|^2)` centered on each cavity grid
!> point, in native cavity order: `<G>`, `<(r-r_i) G>`, `<(r-r_i)(r-r_i) G>` and
!> `<(r-r_i) |r-r_i|^2 G>`. The width `w_i` is the model's own; read the areas
!> back from the live cavity before forming them.
!>
!> They must be rebuilt after every cavity update. Supplying moments whose grid
!> size disagrees with the live cavity is an error rather than a resize: stale
!> moments would otherwise produce a plausible energy for a surface that no
!> longer exists.
   subroutine general_model_supply_gostshyp_api(verror, vmodel, ngrid, c_gt, c_pt, c_mt, c_rt) &
      & bind(C, name=namespace//"general_model_supply_gostshyp")
      !> Error and model handles
      type(c_ptr), value :: verror, vmodel
      !> Cavity grid size
      integer(c_int), value :: ngrid
      !> Input arrays
      type(c_ptr), value :: c_gt, c_pt, c_mt, c_rt
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded model wrapper
      type(vp_model), pointer :: model
      !> Scalar input view
      real(c_double), pointer :: gt(:)
      !> Vector input views
      real(c_double), pointer :: pt(:, :), rt(:, :)
      !> Tensor input view
      real(c_double), pointer :: mt(:, :, :)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmodel) .or. .not. c_associated(c_gt) .or. &
          .not. c_associated(c_pt) .or. .not. c_associated(c_mt) .or. &
          .not. c_associated(c_rt) .or. ngrid <= 0) then
         call api_error(error%ptr, "general_model_supply_gostshyp", &
                        "Invalid model, moment array, or grid size")
         return
      end if
      call c_f_pointer(vmodel, model)
      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "general_model_supply_gostshyp", &
                        "Model is not initialized")
         return
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         if (.not. general%updated .or. ngrid /= general%cavity%ngrid) then
            call api_error(error%ptr, "general_model_supply_gostshyp", &
                           "Grid size does not match the updated general-model cavity")
            return
         end if
      class default
         call api_error(error%ptr, "general_model_supply_gostshyp", &
                        "Model is not a general solvation model")
         return
      end select

      call c_f_pointer(c_gt, gt, [int(ngrid)])
      call c_f_pointer(c_pt, pt, [3, int(ngrid)])
      call c_f_pointer(c_mt, mt, [3, 3, int(ngrid)])
      call c_f_pointer(c_rt, rt, [3, int(ngrid)])
      model%coupling%gostshyp%gt = real(gt, wp)
      model%coupling%gostshyp%pt = real(pt, wp)
      model%coupling%gostshyp%mt = real(mt, wp)
      model%coupling%gostshyp%rt = real(rt, wp)

   end subroutine general_model_supply_gostshyp_api

!> Return the direct trace response from a general model.
!>
!> `surface_charge` may be null, which states that the caller does not want the
!> channel. A non-null pointer is a requirement: if no component produces surface
!> charges the call fails rather than handing back a zero-filled buffer that
!> reads like a converged answer.
   subroutine general_model_get_trace_response_api(verror, vmodel, ngrid, c_surface_charge) &
      & bind(C, name=namespace//"general_model_get_trace_response")
      !> Error and model handles
      type(c_ptr), value :: verror, vmodel
      !> Electrostatic grid size
      integer(c_int), value :: ngrid
      !> Output array
      type(c_ptr), value :: c_surface_charge
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded model wrapper
      type(vp_model), pointer :: model
      !> Output view
      real(c_double), pointer :: surface_charge(:)
      !> Fortran response result
      type(response_type) :: response
      !> Model error
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmodel)) then
         call api_error(error%ptr, "general_model_get_trace_response", "Null pointer provided")
         return
      end if
      call c_f_pointer(vmodel, model)
      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "general_model_get_trace_response", &
                        "Model is not initialized")
         return
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         if (.not. general%updated .or. ngrid /= general%cavity%ngrid) then
            call api_error(error%ptr, "general_model_get_trace_response", &
                           "Grid size does not match the updated general-model cavity")
            return
         end if
      class default
         call api_error(error%ptr, "general_model_get_trace_response", &
                        "Model is not a general solvation model")
         return
      end select

      if (c_associated(c_surface_charge)) then
         call c_f_pointer(c_surface_charge, surface_charge, [int(ngrid)])
         surface_charge = 0.0_c_double
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         call general%get_trace_response(model%coupling, response, model_error)
      class default
         call fatal_error(model_error, "Model is not a general solvation model")
      end select
      if (allocated(model_error)) then
         call api_error(error%ptr, "general_model_get_trace_response", model_error%message)
         return
      end if

      if (c_associated(c_surface_charge)) then
         if (.not. allocated(response%electrostatics%surface_charge)) then
            call api_channel_absent(error%ptr, "general_model_get_trace_response", &
                                    "electrostatics%surface_charge")
            return
         end if
         surface_charge = response%electrostatics%surface_charge
      end if

   end subroutine general_model_get_trace_response_api

!> Return the direct and level-set response channels.
!>
!> Every output pointer is optional. A null pointer states that the caller does
!> not want that channel; a non-null pointer states that the caller requires it,
!> and the call fails if this model configuration produces nothing for it. That
!> distinction matters because "absent" is a legitimate physical answer here --
!> a cavity with field-independent geometry has no level-set response -- and a
!> silently zero-filled buffer is indistinguishable from a converged zero.
   subroutine general_model_get_response_api(verror, vmodel, ngrid, c_surface_charge, &
                                             c_w_value, c_w_gradient, c_w_hessian) &
      & bind(C, name=namespace//"general_model_get_response")
      !> Error and model handles
      type(c_ptr), value :: verror, vmodel
      !> Cavity grid size
      integer(c_int), value :: ngrid
      !> Output arrays, each optional
      type(c_ptr), value :: c_surface_charge, c_w_value, c_w_gradient, c_w_hessian
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded model wrapper
      type(vp_model), pointer :: model
      !> Scalar outputs
      real(c_double), pointer :: surface_charge(:), w_value(:)
      !> Tensor outputs
      real(c_double), pointer :: w_gradient(:, :), w_hessian(:, :, :)
      !> Fortran response result
      type(response_type) :: response
      !> Model error
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmodel)) then
         call api_error(error%ptr, "general_model_get_response", "Null pointer provided")
         return
      end if
      call c_f_pointer(vmodel, model)
      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "general_model_get_response", "Model is not initialized")
         return
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         if (.not. general%updated .or. ngrid /= general%cavity%ngrid) then
            call api_error(error%ptr, "general_model_get_response", &
                           "Grid size does not match the updated general-model cavity")
            return
         end if
      class default
         call api_error(error%ptr, "general_model_get_response", &
                        "Model is not a general solvation model")
         return
      end select

      if (c_associated(c_surface_charge)) then
         call c_f_pointer(c_surface_charge, surface_charge, [int(ngrid)])
         surface_charge = 0.0_c_double
      end if
      if (c_associated(c_w_value)) then
         call c_f_pointer(c_w_value, w_value, [int(ngrid)])
         w_value = 0.0_c_double
      end if
      if (c_associated(c_w_gradient)) then
         call c_f_pointer(c_w_gradient, w_gradient, [3, int(ngrid)])
         w_gradient = 0.0_c_double
      end if
      if (c_associated(c_w_hessian)) then
         call c_f_pointer(c_w_hessian, w_hessian, [3, 3, int(ngrid)])
         w_hessian = 0.0_c_double
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         call general%get_response(model%coupling, response, model_error)
      class default
         call fatal_error(model_error, "Model is not a general solvation model")
      end select
      if (allocated(model_error)) then
         call api_error(error%ptr, "general_model_get_response", model_error%message)
         return
      end if

      if (c_associated(c_surface_charge)) then
         if (.not. allocated(response%electrostatics%surface_charge)) then
            call api_channel_absent(error%ptr, "general_model_get_response", "electrostatics%surface_charge")
            return
         end if
         surface_charge = response%electrostatics%surface_charge
      end if
      if (c_associated(c_w_value)) then
         if (.not. allocated(response%lsf%w_value)) then
            call api_channel_absent(error%ptr, "general_model_get_response", "lsf%w_value")
            return
         end if
         w_value = response%lsf%w_value
      end if
      if (c_associated(c_w_gradient)) then
         if (.not. allocated(response%lsf%w_gradient)) then
            call api_channel_absent(error%ptr, "general_model_get_response", "lsf%w_gradient")
            return
         end if
         w_gradient = response%lsf%w_gradient
      end if
      if (c_associated(c_w_hessian)) then
         if (.not. allocated(response%lsf%w_hessian)) then
            call api_channel_absent(error%ptr, "general_model_get_response", "lsf%w_hessian")
            return
         end if
         w_hessian = response%lsf%w_hessian
      end if

   end subroutine general_model_get_response_api

!> Return every response channel, including the host Gaussian amplitudes.
!>
!> Identical to `general_model_get_response` but also reports the GOSTSHYP
!> amplitudes, which are the component's half of the Fock matrix: the host
!> completes it as
!>
!>    F_uv += sum_i [ w_overlap(i) g_uv,i + w_normal_deriv(i) f_uv,i ]
!>
!> The same optional-pointer contract applies to every channel, so a host
!> without GOSTSHYP passes null for the two amplitude pointers rather than
!> receiving zeros it cannot distinguish from a real result.
!>
!> This exists rather than a standalone amplitude getter because assembling a
!> response runs every component and contracts the cavity surface adjoints
!> once; splitting the read in two would pay that cost twice.
   subroutine general_model_get_response_extended_api(verror, vmodel, ngrid, c_surface_charge, &
                                                      c_w_value, c_w_gradient, c_w_hessian, &
                                                      c_w_overlap, c_w_normal_deriv) &
      & bind(C, name=namespace//"general_model_get_response_extended")
      !> Error and model handles
      type(c_ptr), value :: verror, vmodel
      !> Cavity grid size
      integer(c_int), value :: ngrid
      !> Output arrays, each optional
      type(c_ptr), value :: c_surface_charge, c_w_value, c_w_gradient, c_w_hessian
      !> Optional Gaussian-amplitude outputs
      type(c_ptr), value :: c_w_overlap, c_w_normal_deriv
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded model wrapper
      type(vp_model), pointer :: model
      !> Scalar outputs
      real(c_double), pointer :: surface_charge(:), w_value(:)
      !> Tensor outputs
      real(c_double), pointer :: w_gradient(:, :), w_hessian(:, :, :)
      !> Gaussian-amplitude outputs
      real(c_double), pointer :: w_overlap(:), w_normal_deriv(:)
      !> Fortran response result
      type(response_type) :: response
      !> Model error
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmodel)) then
         call api_error(error%ptr, "general_model_get_response_extended", "Null pointer provided")
         return
      end if
      call c_f_pointer(vmodel, model)
      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "general_model_get_response_extended", "Model is not initialized")
         return
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         if (.not. general%updated .or. ngrid /= general%cavity%ngrid) then
            call api_error(error%ptr, "general_model_get_response_extended", &
                           "Grid size does not match the updated general-model cavity")
            return
         end if
      class default
         call api_error(error%ptr, "general_model_get_response_extended", &
                        "Model is not a general solvation model")
         return
      end select

      if (c_associated(c_surface_charge)) then
         call c_f_pointer(c_surface_charge, surface_charge, [int(ngrid)])
         surface_charge = 0.0_c_double
      end if
      if (c_associated(c_w_value)) then
         call c_f_pointer(c_w_value, w_value, [int(ngrid)])
         w_value = 0.0_c_double
      end if
      if (c_associated(c_w_gradient)) then
         call c_f_pointer(c_w_gradient, w_gradient, [3, int(ngrid)])
         w_gradient = 0.0_c_double
      end if
      if (c_associated(c_w_hessian)) then
         call c_f_pointer(c_w_hessian, w_hessian, [3, 3, int(ngrid)])
         w_hessian = 0.0_c_double
      end if
      if (c_associated(c_w_overlap)) then
         call c_f_pointer(c_w_overlap, w_overlap, [int(ngrid)])
         w_overlap = 0.0_c_double
      end if
      if (c_associated(c_w_normal_deriv)) then
         call c_f_pointer(c_w_normal_deriv, w_normal_deriv, [int(ngrid)])
         w_normal_deriv = 0.0_c_double
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         call general%get_response(model%coupling, response, model_error)
      class default
         call fatal_error(model_error, "Model is not a general solvation model")
      end select
      if (allocated(model_error)) then
         call api_error(error%ptr, "general_model_get_response_extended", model_error%message)
         return
      end if

      if (c_associated(c_surface_charge)) then
         if (.not. allocated(response%electrostatics%surface_charge)) then
            call api_channel_absent(error%ptr, "general_model_get_response_extended", "electrostatics%surface_charge")
            return
         end if
         surface_charge = response%electrostatics%surface_charge
      end if
      if (c_associated(c_w_value)) then
         if (.not. allocated(response%lsf%w_value)) then
            call api_channel_absent(error%ptr, "general_model_get_response_extended", "lsf%w_value")
            return
         end if
         w_value = response%lsf%w_value
      end if
      if (c_associated(c_w_gradient)) then
         if (.not. allocated(response%lsf%w_gradient)) then
            call api_channel_absent(error%ptr, "general_model_get_response_extended", "lsf%w_gradient")
            return
         end if
         w_gradient = response%lsf%w_gradient
      end if
      if (c_associated(c_w_hessian)) then
         if (.not. allocated(response%lsf%w_hessian)) then
            call api_channel_absent(error%ptr, "general_model_get_response_extended", "lsf%w_hessian")
            return
         end if
         w_hessian = response%lsf%w_hessian
      end if
      if (c_associated(c_w_overlap)) then
         if (.not. allocated(response%gostshyp%w_overlap)) then
            call api_channel_absent(error%ptr, "general_model_get_response_extended", "gostshyp%w_overlap")
            return
         end if
         w_overlap = response%gostshyp%w_overlap
      end if
      if (c_associated(c_w_normal_deriv)) then
         if (.not. allocated(response%gostshyp%w_normal_deriv)) then
            call api_channel_absent(error%ptr, "general_model_get_response_extended", "gostshyp%w_normal_deriv")
            return
         end if
         w_normal_deriv = response%gostshyp%w_normal_deriv
      end if

   end subroutine general_model_get_response_extended_api

!> Return the general-model nuclear gradient.
   subroutine general_model_get_gradient_api(verror, vmodel, nat, c_gradient) &
         & bind(C, name=namespace//"general_model_get_gradient")
      !> Error and model handles
      type(c_ptr), value :: verror, vmodel
      !> Number of atoms
      integer(c_int), value :: nat
      !> Output gradient
      type(c_ptr), value :: c_gradient
      !> Decoded error wrapper
      type(vp_error), pointer :: error
      !> Decoded model wrapper
      type(vp_model), pointer :: model
      !> Output view
      real(c_double), pointer :: gradient(:, :)
      !> Fortran gradient
      real(wp), allocatable :: local(:, :)
      !> Model error
      type(error_type), allocatable :: model_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vmodel) .or. .not. c_associated(c_gradient)) then
         call api_error(error%ptr, "general_model_get_gradient", "Null pointer provided")
         return
      end if
      call c_f_pointer(vmodel, model)
      if (.not. allocated(model%ptr)) then
         call api_error(error%ptr, "general_model_get_gradient", "Model is not initialized")
         return
      end if

      select type (general => model%ptr)
      type is (solvation_model_general)
         if (.not. general%updated .or. nat /= general%cavity%nsph) then
            call api_error(error%ptr, "general_model_get_gradient", &
                           "Atom count does not match the updated general-model cavity")
            return
         end if
      class default
         call api_error(error%ptr, "general_model_get_gradient", &
                        "Model is not a general solvation model")
         return
      end select

      call c_f_pointer(c_gradient, gradient, [3, int(nat)])
      allocate (local(3, int(nat)), source=0.0_wp)
      select type (general => model%ptr)
      type is (solvation_model_general)
         call general%get_gradient(model%coupling, local, model_error)
      class default
         call fatal_error(model_error, "Model is not a general solvation model")
      end select
      if (allocated(model_error)) then
         call api_error(error%ptr, "general_model_get_gradient", model_error%message)
         return
      end if
      gradient = real(local, c_double)

   end subroutine general_model_get_gradient_api

!> Decode the optional master tolerance handed to a DROP cavity constructor.
!> A null pointer selects the compiled DROP default, so every entry point can
!> forward one unconditional `tolerance=` argument to `new_cavity_drop` instead
!> of branching on presence at each construction site.
!> @param[in]  c_tolerance Caller pointer to the master tolerance, or NULL
!> @param[out] tolerance   Decoded tolerance (compiled default when NULL)
!> @return .true. when the decoded value is usable, .false. if non-positive
   logical function decode_drop_tolerance(c_tolerance, tolerance) result(ok)
      !> Caller-supplied tolerance pointer (may be NULL)
      type(c_ptr), intent(in) :: c_tolerance
      !> Decoded master tolerance
      real(wp), intent(out) :: tolerance
      !> Decoded caller value
      real(c_double), pointer :: p_tolerance
      !> Default-initialized DROP parameters; single source of the API default
      type(moist_cavity_drop_parameters_type) :: drop_defaults

      tolerance = drop_defaults%tolerance
      ok = .true.
      if (.not. c_associated(c_tolerance)) return

      call c_f_pointer(c_tolerance, p_tolerance)
      tolerance = real(p_tolerance, wp)
      ok = tolerance > 0.0_wp
   end function decode_drop_tolerance

!> Decode the optional controls shared by every DROP constructor.
   subroutine decode_drop_controls(c_proj_maxiter, c_proj_level, &
                                   c_branch_weight_s, c_rho_grid_h, &
                                   c_wleb_prune_level, use_proj_maxiter, &
                                   use_proj_level, use_branch_weight_s, &
                                   use_rho_grid_h, use_wleb_prune_level)
      type(c_ptr), intent(in) :: c_proj_maxiter
      type(c_ptr), intent(in) :: c_proj_level
      type(c_ptr), intent(in) :: c_branch_weight_s
      type(c_ptr), intent(in) :: c_rho_grid_h
      type(c_ptr), intent(in) :: c_wleb_prune_level
      integer, intent(out) :: use_proj_maxiter
      integer, intent(out) :: use_proj_level
      real(wp), intent(out) :: use_branch_weight_s
      real(wp), intent(out) :: use_rho_grid_h
      integer, intent(out) :: use_wleb_prune_level
      integer(c_int), pointer :: p_proj_maxiter, p_proj_level
      integer(c_int), pointer :: p_wleb_prune_level
      real(c_double), pointer :: p_branch_weight_s, p_rho_grid_h
      type(moist_cavity_drop_parameters_type) :: defaults

      use_proj_maxiter = defaults%proj_maxiter
      if (c_associated(c_proj_maxiter)) then
         call c_f_pointer(c_proj_maxiter, p_proj_maxiter)
         use_proj_maxiter = p_proj_maxiter
      end if
      use_proj_level = defaults%proj_level
      if (c_associated(c_proj_level)) then
         call c_f_pointer(c_proj_level, p_proj_level)
         use_proj_level = p_proj_level
      end if
      ! new_cavity_drop has historically used 0.05 as its public default.
      use_branch_weight_s = 0.05_wp
      if (c_associated(c_branch_weight_s)) then
         call c_f_pointer(c_branch_weight_s, p_branch_weight_s)
         use_branch_weight_s = p_branch_weight_s
      end if
      use_rho_grid_h = defaults%rho_grid_h
      if (c_associated(c_rho_grid_h)) then
         call c_f_pointer(c_rho_grid_h, p_rho_grid_h)
         use_rho_grid_h = p_rho_grid_h
      end if
      use_wleb_prune_level = defaults%wleb_prune_level
      if (c_associated(c_wleb_prune_level)) then
         call c_f_pointer(c_wleb_prune_level, p_wleb_prune_level)
         use_wleb_prune_level = p_wleb_prune_level
      end if
   end subroutine decode_drop_controls

!> Allocate and configure a DROP cavity from an already configured LSF.
!> This is the common ownership/error-handling seam for the SvdW and CFC
!> constructors; each public constructor only decodes its own shape settings.
   subroutine new_drop_cavity_from_lsf(error, vradii, nleb, use_debug, use_verbose, &
                                       use_do_fine, use_tolerance, use_proj_maxiter, &
                                       use_proj_level, use_branch_weight_s, use_rho_grid_h, &
                                       use_wleb_prune_level, lsf_template, &
                                       routine_name, vcav)
      type(vp_error), intent(inout) :: error
      type(c_ptr), intent(in) :: vradii
      type(c_ptr), intent(in) :: nleb
      logical, intent(in) :: use_debug
      integer, intent(in) :: use_verbose
      logical, intent(in) :: use_do_fine
      real(wp), intent(in) :: use_tolerance
      integer, intent(in) :: use_proj_maxiter, use_proj_level
      real(wp), intent(in) :: use_branch_weight_s, use_rho_grid_h
      integer, intent(in) :: use_wleb_prune_level
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_template
      character(len=*), intent(in) :: routine_name
      type(c_ptr), intent(out) :: vcav
      type(vp_radii), pointer :: radii_handle
      integer(c_int), pointer :: pnleb
      type(vp_cavity), pointer :: cav
      logical :: use_explicit_radii
      type(error_type), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model

      vcav = c_null_ptr
      use_explicit_radii = c_associated(vradii)
      radii_handle => null()
      if (use_explicit_radii) then
         call c_f_pointer(vradii, radii_handle)
         if (.not. allocated(radii_handle%ptr)) then
            call api_error(error%ptr, routine_name, "Radii handle is not initialized")
            return
         end if
      else
         call new_radii("cpcm", radius_model, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, routine_name, cavity_error%message)
            return
         end if
      end if

      allocate (cav)
      allocate (cavity_type_drop :: cav%ptr)
      call new_context(cav%ctx, verbosity=use_verbose, debug=use_debug)
      pnleb => null()
      if (c_associated(nleb)) call c_f_pointer(nleb, pnleb)

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         if (associated(pnleb)) then
            if (use_explicit_radii) then
               call new_cavity_drop(cavity, cav%ctx, nleb=pnleb, &
                                    tolerance=use_tolerance, do_fine=use_do_fine, &
                                    proj_maxiter=use_proj_maxiter, proj_level=use_proj_level, &
                                    branch_weight_s=use_branch_weight_s, &
                                    rho_grid_h=use_rho_grid_h, &
                                    wleb_prune_level=use_wleb_prune_level, &
                                    radius_model=radii_handle%ptr, lsf_model=lsf_template, &
                                    error=cavity_error)
            else
               call new_cavity_drop(cavity, cav%ctx, nleb=pnleb, &
                                    tolerance=use_tolerance, do_fine=use_do_fine, &
                                    proj_maxiter=use_proj_maxiter, proj_level=use_proj_level, &
                                    branch_weight_s=use_branch_weight_s, &
                                    rho_grid_h=use_rho_grid_h, &
                                    wleb_prune_level=use_wleb_prune_level, &
                                    radius_model=radius_model, lsf_model=lsf_template, &
                                    error=cavity_error)
            end if
         else
            if (use_explicit_radii) then
               call new_cavity_drop(cavity, cav%ctx, tolerance=use_tolerance, &
                                    do_fine=use_do_fine, radius_model=radii_handle%ptr, &
                                    proj_maxiter=use_proj_maxiter, proj_level=use_proj_level, &
                                    branch_weight_s=use_branch_weight_s, &
                                    rho_grid_h=use_rho_grid_h, &
                                    wleb_prune_level=use_wleb_prune_level, &
                                    lsf_model=lsf_template, error=cavity_error)
            else
               call new_cavity_drop(cavity, cav%ctx, tolerance=use_tolerance, &
                                    do_fine=use_do_fine, radius_model=radius_model, &
                                    proj_maxiter=use_proj_maxiter, proj_level=use_proj_level, &
                                    branch_weight_s=use_branch_weight_s, &
                                    rho_grid_h=use_rho_grid_h, &
                                    wleb_prune_level=use_wleb_prune_level, &
                                    lsf_model=lsf_template, error=cavity_error)
            end if
         end if
      end select
      if (allocated(cavity_error)) then
         call api_error(error%ptr, routine_name, cavity_error%message)
         if (associated(cav%ptr)) deallocate (cav%ptr)
         call cav%ctx%delete()
         deallocate (cav)
         return
      end if

      vcav = c_loc(cav)
   end subroutine new_drop_cavity_from_lsf

!> Internal helper to create DROP cavity handles.
   subroutine new_drop_cavity_common(verror, vradii, nleb, c_debug, c_verbose, &
                                     c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine, &
                                     c_tolerance, c_proj_maxiter, c_proj_level, &
                                     c_branch_weight_s, c_rho_grid_h, &
                                     c_wleb_prune_level, routine_name, vcav)
      type(c_ptr), value :: verror
      type(c_ptr), value :: vradii
      type(c_ptr), value :: nleb
      type(c_ptr), value :: c_debug
      type(c_ptr), value :: c_verbose
      type(c_ptr), value :: c_blendk
      type(c_ptr), value :: c_blend1b
      type(c_ptr), value :: c_blend2b
      type(c_ptr), value :: c_blend3b
      type(c_ptr), value :: c_do_fine
      !> Optional master numerical tolerance (NULL selects the DROP default)
      type(c_ptr), value :: c_tolerance
      type(c_ptr), value :: c_proj_maxiter
      type(c_ptr), value :: c_proj_level
      type(c_ptr), value :: c_branch_weight_s
      type(c_ptr), value :: c_rho_grid_h
      type(c_ptr), value :: c_wleb_prune_level
      character(len=*), intent(in) :: routine_name
      type(c_ptr), intent(out) :: vcav
      type(vp_error), pointer :: error
      logical(c_bool), pointer :: p_debug
      integer(c_int), pointer :: p_verbose
      real(c_double), pointer :: p_blendk
      real(c_double), pointer :: p_blend1b
      real(c_double), pointer :: p_blend2b
      real(c_double), pointer :: p_blend3b
      logical(c_bool), pointer :: p_do_fine
      logical :: use_debug
      integer :: use_verbose
      real(wp) :: use_blendk, use_blend1b, use_blend2b, use_blend3b
      !> Default-initialized SvdW shape parameters.
      !>
      !> The fallbacks below are read from this rather than written out again, so
      !> a caller passing NULL through the C API lands on exactly the surface a
      !> Fortran caller gets by not overriding. Duplicating the literals here is
      !> how the two entry points silently drifted apart before.
      type(moist_cavity_drop_lsf_svdw_param_type) :: svdw_defaults
      logical :: use_do_fine
      real(wp) :: use_tolerance
      integer :: use_proj_maxiter, use_proj_level, use_wleb_prune_level
      real(wp) :: use_branch_weight_s, use_rho_grid_h
      type(moist_cavity_drop_lsf_svdw_type) :: svdw_template

      vcav = c_null_ptr

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      ! Parse optional debug flag (default: false)
      p_debug => null()
      if (c_associated(c_debug)) then
         call c_f_pointer(c_debug, p_debug)
         use_debug = p_debug
      else
         use_debug = .false.
      end if

      ! Parse optional verbose flag (default: 0)
      p_verbose => null()
      if (c_associated(c_verbose)) then
         call c_f_pointer(c_verbose, p_verbose)
         use_verbose = p_verbose
      else
         use_verbose = 0
      end if

      ! Parse optional SvdW shape parameters from their native defaults.
      p_blendk => null()
      if (c_associated(c_blendk)) then
         call c_f_pointer(c_blendk, p_blendk)
         use_blendk = p_blendk
      else
         use_blendk = svdw_defaults%blend_k
      end if

      p_blend1b => null()
      if (c_associated(c_blend1b)) then
         call c_f_pointer(c_blend1b, p_blend1b)
         use_blend1b = p_blend1b
      else
         use_blend1b = svdw_defaults%blend_1b
      end if

      p_blend2b => null()
      if (c_associated(c_blend2b)) then
         call c_f_pointer(c_blend2b, p_blend2b)
         use_blend2b = p_blend2b
      else
         use_blend2b = svdw_defaults%blend_2b
      end if

      p_blend3b => null()
      if (c_associated(c_blend3b)) then
         call c_f_pointer(c_blend3b, p_blend3b)
         use_blend3b = p_blend3b
      else
         use_blend3b = svdw_defaults%blend_3b
      end if

      ! Parse optional do_fine flag (default: false)
      p_do_fine => null()
      if (c_associated(c_do_fine)) then
         call c_f_pointer(c_do_fine, p_do_fine)
         use_do_fine = p_do_fine
      else
         use_do_fine = .false.
      end if

      ! Parse optional master tolerance (default: compiled DROP default)
      if (.not. decode_drop_tolerance(c_tolerance, use_tolerance)) then
         call api_error(error%ptr, routine_name, "DROP tolerance must be positive")
         return
      end if
      call decode_drop_controls(c_proj_maxiter, c_proj_level, &
                                c_branch_weight_s, c_rho_grid_h, c_wleb_prune_level, &
                                use_proj_maxiter, use_proj_level, use_branch_weight_s, &
                                use_rho_grid_h, use_wleb_prune_level)

      ! The cavity derives the LSF screening threshold from its own tolerance,
      ! so only the shape parameters are configured here.
      call svdw_template%new(blend_k=use_blendk, blend_1b=use_blend1b, &
                             blend_2b=use_blend2b, blend_3b=use_blend3b)
      call new_drop_cavity_from_lsf(error, vradii, nleb, use_debug, use_verbose, &
                                    use_do_fine, use_tolerance, use_proj_maxiter, &
                                    use_proj_level, use_branch_weight_s, use_rho_grid_h, &
                                    use_wleb_prune_level, svdw_template, &
                                    routine_name, vcav)
   end subroutine new_drop_cavity_common

!> Create new DROP cavity handle with default CPCM radii model.
   function new_drop_cavity_api(verror, nleb, c_debug, c_verbose, &
         c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine, c_tolerance, &
         c_proj_maxiter, c_proj_level, c_branch_weight_s, c_rho_grid_h, &
         c_wleb_prune_level) result(vcav) &
         & bind(C, name=namespace//"new_drop_cavity")
      type(c_ptr), value :: verror
      type(c_ptr), value :: nleb
      type(c_ptr), value :: c_debug
      type(c_ptr), value :: c_verbose
      type(c_ptr), value :: c_blendk
      type(c_ptr), value :: c_blend1b
      type(c_ptr), value :: c_blend2b
      type(c_ptr), value :: c_blend3b
      type(c_ptr), value :: c_do_fine
      !> Optional master numerical tolerance (NULL selects the DROP default)
      type(c_ptr), value :: c_tolerance
      type(c_ptr), value :: c_proj_maxiter
      type(c_ptr), value :: c_proj_level
      type(c_ptr), value :: c_branch_weight_s
      type(c_ptr), value :: c_rho_grid_h
      type(c_ptr), value :: c_wleb_prune_level
      type(c_ptr) :: vcav

      call new_drop_cavity_common(verror, c_null_ptr, nleb, c_debug, c_verbose, &
                                  c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine, &
                                  c_tolerance, c_proj_maxiter, c_proj_level, &
                                  c_branch_weight_s, c_rho_grid_h, c_wleb_prune_level, &
                                  "new_drop_cavity_api", vcav)
   end function new_drop_cavity_api

!> Create new DROP cavity handle with an explicit radii model.
   function new_drop_cavity_with_radii_api(verror, vradii, nleb, c_debug, c_verbose, &
         c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine, c_tolerance, &
         c_proj_maxiter, c_proj_level, c_branch_weight_s, c_rho_grid_h, &
         c_wleb_prune_level) result(vcav) &
         & bind(C, name=namespace//"new_drop_cavity_with_radii")
      type(c_ptr), value :: verror
      type(c_ptr), value :: vradii
      type(c_ptr), value :: nleb
      type(c_ptr), value :: c_debug
      type(c_ptr), value :: c_verbose
      type(c_ptr), value :: c_blendk
      type(c_ptr), value :: c_blend1b
      type(c_ptr), value :: c_blend2b
      type(c_ptr), value :: c_blend3b
      type(c_ptr), value :: c_do_fine
      !> Optional master numerical tolerance (NULL selects the DROP default)
      type(c_ptr), value :: c_tolerance
      type(c_ptr), value :: c_proj_maxiter
      type(c_ptr), value :: c_proj_level
      type(c_ptr), value :: c_branch_weight_s
      type(c_ptr), value :: c_rho_grid_h
      type(c_ptr), value :: c_wleb_prune_level
      type(c_ptr) :: vcav

      call new_drop_cavity_common(verror, vradii, nleb, c_debug, c_verbose, &
                                  c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine, &
                                  c_tolerance, c_proj_maxiter, c_proj_level, &
                                  c_branch_weight_s, c_rho_grid_h, c_wleb_prune_level, &
                                  "new_drop_cavity_with_radii_api", vcav)
   end function new_drop_cavity_with_radii_api

!> Create a CFC-DROP cavity handle with default CPCM radii.
   function new_cfc_drop_cavity_api(verror, nleb, c_debug, c_verbose, &
         c_a1, c_a2, c_c, c_m, c_screen_k, c_do_fine, c_tolerance, &
         c_proj_maxiter, c_proj_level, c_branch_weight_s, c_rho_grid_h, &
         c_wleb_prune_level) result(vcav) &
         & bind(C, name=namespace//"new_cfc_drop_cavity")
      type(c_ptr), value :: verror
      type(c_ptr), value :: nleb
      type(c_ptr), value :: c_debug
      type(c_ptr), value :: c_verbose
      type(c_ptr), value :: c_a1
      type(c_ptr), value :: c_a2
      type(c_ptr), value :: c_c
      type(c_ptr), value :: c_m
      type(c_ptr), value :: c_screen_k
      type(c_ptr), value :: c_do_fine
      type(c_ptr), value :: c_tolerance
      type(c_ptr), value :: c_proj_maxiter
      type(c_ptr), value :: c_proj_level
      type(c_ptr), value :: c_branch_weight_s
      type(c_ptr), value :: c_rho_grid_h
      type(c_ptr), value :: c_wleb_prune_level
      type(c_ptr) :: vcav
      type(vp_error), pointer :: error
      logical(c_bool), pointer :: p_debug, p_do_fine
      integer(c_int), pointer :: p_verbose, p_m
      real(c_double), pointer :: p_a1, p_a2, p_c, p_screen_k
      logical :: use_debug, use_do_fine
      integer :: use_verbose, use_m
      real(wp) :: use_a1, use_a2, use_c, use_screen_k, use_tolerance
      integer :: use_proj_maxiter, use_proj_level, use_wleb_prune_level
      real(wp) :: use_branch_weight_s, use_rho_grid_h
      type(moist_cavity_drop_lsf_cfc_type) :: cfc_template

      vcav = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      use_debug = .false.
      if (c_associated(c_debug)) then
         call c_f_pointer(c_debug, p_debug)
         use_debug = p_debug
      end if
      use_verbose = 0
      if (c_associated(c_verbose)) then
         call c_f_pointer(c_verbose, p_verbose)
         use_verbose = p_verbose
      end if
      use_do_fine = .false.
      if (c_associated(c_do_fine)) then
         call c_f_pointer(c_do_fine, p_do_fine)
         use_do_fine = p_do_fine
      end if

      use_a1 = cfc_template%a1
      if (c_associated(c_a1)) then
         call c_f_pointer(c_a1, p_a1)
         use_a1 = p_a1
      end if
      use_a2 = cfc_template%a2
      if (c_associated(c_a2)) then
         call c_f_pointer(c_a2, p_a2)
         use_a2 = p_a2
      end if
      use_c = cfc_template%c
      if (c_associated(c_c)) then
         call c_f_pointer(c_c, p_c)
         use_c = p_c
      end if
      use_m = cfc_template%m
      if (c_associated(c_m)) then
         call c_f_pointer(c_m, p_m)
         use_m = p_m
      end if
      use_screen_k = cfc_template%screen_k
      if (c_associated(c_screen_k)) then
         call c_f_pointer(c_screen_k, p_screen_k)
         use_screen_k = p_screen_k
      end if
      if (.not. decode_drop_tolerance(c_tolerance, use_tolerance)) then
         call api_error(error%ptr, "new_cfc_drop_cavity_api", &
                        "DROP tolerance must be positive")
         return
      end if
      call decode_drop_controls(c_proj_maxiter, c_proj_level, &
                                c_branch_weight_s, c_rho_grid_h, c_wleb_prune_level, &
                                use_proj_maxiter, use_proj_level, use_branch_weight_s, &
                                use_rho_grid_h, use_wleb_prune_level)

      call cfc_template%new(a1=use_a1, a2=use_a2, c=use_c, m=use_m, &
                            screen_k=use_screen_k)
      call new_drop_cavity_from_lsf(error, c_null_ptr, nleb, use_debug, use_verbose, &
                                    use_do_fine, use_tolerance, use_proj_maxiter, &
                                    use_proj_level, use_branch_weight_s, use_rho_grid_h, &
                                    use_wleb_prune_level, cfc_template, &
                                    "new_cfc_drop_cavity_api", vcav)
   end function new_cfc_drop_cavity_api

!> Create an iSwiG cavity handle with default CPCM radii.
   function new_iswig_cavity_api(verror, nleb, c_debug, c_verbose, &
         c_cut_a, c_cut_f) result(vcav) &
         & bind(C, name=namespace//"new_iswig_cavity")
      type(c_ptr), value :: verror
      type(c_ptr), value :: nleb
      type(c_ptr), value :: c_debug
      type(c_ptr), value :: c_verbose
      type(c_ptr), value :: c_cut_a
      type(c_ptr), value :: c_cut_f
      type(c_ptr) :: vcav
      type(vp_error), pointer :: error
      integer(c_int), pointer :: p_nleb, p_verbose
      logical(c_bool), pointer :: p_debug
      real(c_double), pointer :: p_cut_a, p_cut_f
      integer :: use_nleb, use_verbose
      logical :: use_debug
      real(wp) :: use_cut_a, use_cut_f
      type(cavity_type_iswig) :: iswig_defaults
      type(vp_cavity), pointer :: cav
      class(radius_type), allocatable :: radius_model
      type(error_type), allocatable :: cavity_error

      vcav = c_null_ptr
      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      use_nleb = iswig_defaults%num_leb
      if (c_associated(nleb)) then
         call c_f_pointer(nleb, p_nleb)
         use_nleb = p_nleb
      end if
      use_debug = .false.
      if (c_associated(c_debug)) then
         call c_f_pointer(c_debug, p_debug)
         use_debug = p_debug
      end if
      use_verbose = 0
      if (c_associated(c_verbose)) then
         call c_f_pointer(c_verbose, p_verbose)
         use_verbose = p_verbose
      end if
      use_cut_a = iswig_defaults%cut_a
      if (c_associated(c_cut_a)) then
         call c_f_pointer(c_cut_a, p_cut_a)
         use_cut_a = p_cut_a
      end if
      use_cut_f = iswig_defaults%cut_f
      if (c_associated(c_cut_f)) then
         call c_f_pointer(c_cut_f, p_cut_f)
         use_cut_f = p_cut_f
      end if

      call new_radii("cpcm", radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call api_error(error%ptr, "new_iswig_cavity_api", cavity_error%message)
         return
      end if
      allocate (cav)
      allocate (cavity_type_iswig :: cav%ptr)
      call new_context(cav%ctx, verbosity=use_verbose, debug=use_debug)
      select type (cavity => cav%ptr)
      type is (cavity_type_iswig)
         call new_cavity_iswig(cavity, cav%ctx, nleb=use_nleb, &
                               cut_a=use_cut_a, cut_f=use_cut_f, &
                               radius_model=radius_model, error=cavity_error)
      end select
      if (allocated(cavity_error)) then
         call api_error(error%ptr, "new_iswig_cavity_api", cavity_error%message)
         if (associated(cav%ptr)) deallocate (cav%ptr)
         call cav%ctx%delete()
         deallocate (cav)
         return
      end if

      vcav = c_loc(cav)
   end function new_iswig_cavity_api

!> Create new DROP cavity handle using an external isodensity LSF callback.
!>
!> The callback must follow the `moist_isodensity_lsf_callback` ABI: return 0
!> once every requested buffer is written, or any nonzero status to report that
!> it could not evaluate. A nonzero status aborts the build and surfaces from
!> `update_cavity` as an ordinary API error naming that status; no cavity data
!> are produced.
   function new_drop_cavity_isodensity_callback_api(verror, callback, context, &
         c_scale, nleb, c_debug, c_verbose, c_do_fine, c_wleb_prune_level, &
         c_tolerance) result(vcav) &
         & bind(C, name=namespace//"new_drop_cavity_isodensity_callback")
      type(c_ptr), value :: verror
      type(c_funptr), value :: callback
      type(c_ptr), value :: context
      type(c_ptr), value :: c_scale
      type(c_ptr), value :: nleb
      type(c_ptr), value :: c_debug
      type(c_ptr), value :: c_verbose
      type(c_ptr), value :: c_do_fine
      type(c_ptr), value :: c_wleb_prune_level
      !> Optional master numerical tolerance (NULL selects the DROP default)
      type(c_ptr), value :: c_tolerance
      type(c_ptr) :: vcav

      type(vp_error), pointer :: error
      integer(c_int), pointer :: pnleb
      real(c_double), pointer :: p_scale
      logical(c_bool), pointer :: p_debug
      integer(c_int), pointer :: p_verbose
      logical(c_bool), pointer :: p_do_fine
      integer(c_int), pointer :: p_wleb_prune_level
      type(vp_cavity), pointer :: cav
      real(wp) :: use_scale
      logical :: use_debug
      integer :: use_verbose
      logical :: use_do_fine
      integer :: use_wleb_prune_level
      real(wp) :: use_tolerance
      type(error_type), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model

      vcav = c_null_ptr

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(callback)) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_callback_api", &
                        "Isodensity callback is missing")
         return
      end if

      p_scale => null()
      if (c_associated(c_scale)) then
         call c_f_pointer(c_scale, p_scale)
         use_scale = p_scale
      else
         use_scale = 1000.0_wp
      end if

      p_debug => null()
      if (c_associated(c_debug)) then
         call c_f_pointer(c_debug, p_debug)
         use_debug = p_debug
      else
         use_debug = .false.
      end if

      p_verbose => null()
      if (c_associated(c_verbose)) then
         call c_f_pointer(c_verbose, p_verbose)
         use_verbose = p_verbose
      else
         use_verbose = 0
      end if

      p_do_fine => null()
      if (c_associated(c_do_fine)) then
         call c_f_pointer(c_do_fine, p_do_fine)
         use_do_fine = p_do_fine
      else
         use_do_fine = .false.
      end if

      p_wleb_prune_level => null()
      if (c_associated(c_wleb_prune_level)) then
         call c_f_pointer(c_wleb_prune_level, p_wleb_prune_level)
         use_wleb_prune_level = p_wleb_prune_level
      else
         use_wleb_prune_level = 0
      end if

      call new_radii("cpcm", radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_callback_api", &
                        cavity_error%message)
         return
      end if

      ! Parse optional master tolerance (default: compiled DROP default)
      if (.not. decode_drop_tolerance(c_tolerance, use_tolerance)) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_callback_api", &
                        "DROP tolerance must be positive")
         return
      end if

      allocate (cav)
      allocate (cavity_type_drop :: cav%ptr)
      call new_context(cav%ctx, verbosity=use_verbose, debug=use_debug)
      pnleb => null()
      if (c_associated(nleb)) call c_f_pointer(nleb, pnleb)

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         block
            type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf_template

            call lsf_template%new(callback, context, scale=use_scale)

            if (associated(pnleb)) then
               call new_cavity_drop(cavity, cav%ctx, nleb=pnleb, &
                                    tolerance=use_tolerance, &
                                    do_fine=use_do_fine, &
                                    wleb_prune_level=use_wleb_prune_level, &
                                    radius_model=radius_model, lsf_model=lsf_template, &
                                    error=cavity_error)
            else
               call new_cavity_drop(cavity, cav%ctx, &
                                    tolerance=use_tolerance, &
                                    do_fine=use_do_fine, wleb_prune_level=use_wleb_prune_level, &
                                    radius_model=radius_model, &
                                    lsf_model=lsf_template, error=cavity_error)
            end if
         end block
      end select

      if (allocated(cavity_error)) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_callback_api", &
                        cavity_error%message)
         if (associated(cav%ptr)) deallocate (cav%ptr)
         deallocate (cav)
         return
      end if

      vcav = c_loc(cav)
   end function new_drop_cavity_isodensity_callback_api

!> Create new DROP cavity handle using the internal isodensity LSF.
!>
!> The cartesian-monomial Gaussian basis is supplied here; the (transformed)
!> density matrix is installed separately via set_isodensity_density before each
!> cavity build.  The cartesian-component ordering can be queried through
!> get_isodensity_cart_layout so the host can build a matching density transform.
   function new_drop_cavity_isodensity_internal_api(verror, nshell, c_shell_atom, &
         c_shell_l, c_shell_nprim, c_exps, c_coeffs, rho_iso, c_scale, nleb, &
         c_debug, c_verbose, c_do_fine, c_wleb_prune_level, c_tolerance) result(vcav) &
         & bind(C, name=namespace//"new_drop_cavity_isodensity_internal")
      type(c_ptr), value :: verror
      integer(c_int), value :: nshell
      type(c_ptr), value :: c_shell_atom
      type(c_ptr), value :: c_shell_l
      type(c_ptr), value :: c_shell_nprim
      type(c_ptr), value :: c_exps
      type(c_ptr), value :: c_coeffs
      real(c_double), value :: rho_iso
      type(c_ptr), value :: c_scale
      type(c_ptr), value :: nleb
      type(c_ptr), value :: c_debug
      type(c_ptr), value :: c_verbose
      type(c_ptr), value :: c_do_fine
      type(c_ptr), value :: c_wleb_prune_level
      !> Optional master numerical tolerance (NULL selects the DROP default)
      type(c_ptr), value :: c_tolerance
      type(c_ptr) :: vcav

      type(vp_error), pointer :: error
      integer(c_int), pointer :: shell_atom(:), shell_l(:), shell_nprim(:)
      real(c_double), pointer :: exps(:), coeffs(:)
      integer(c_int), pointer :: pnleb
      real(c_double), pointer :: p_scale
      logical(c_bool), pointer :: p_debug
      integer(c_int), pointer :: p_verbose
      logical(c_bool), pointer :: p_do_fine
      integer(c_int), pointer :: p_wleb_prune_level
      type(vp_cavity), pointer :: cav
      real(wp) :: use_scale
      logical :: use_debug
      integer :: use_verbose
      logical :: use_do_fine
      integer :: use_wleb_prune_level
      integer, allocatable :: sh_atom(:), sh_l(:), sh_nprim(:)
      integer :: nprim_tot
      real(wp) :: use_tolerance
      type(error_type), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model

      vcav = c_null_ptr

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (nshell < 1) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_internal_api", &
                        "Number of shells must be positive")
         return
      end if
      if (.not. c_associated(c_shell_atom) .or. .not. c_associated(c_shell_l) &
          .or. .not. c_associated(c_shell_nprim) .or. .not. c_associated(c_exps) &
          .or. .not. c_associated(c_coeffs)) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_internal_api", &
                        "Missing basis array pointer")
         return
      end if

      call c_f_pointer(c_shell_atom, shell_atom, [int(nshell)])
      call c_f_pointer(c_shell_l, shell_l, [int(nshell)])
      call c_f_pointer(c_shell_nprim, shell_nprim, [int(nshell)])

      ! Convert the host's 0-based atom indices to moist's 1-based convention.
      sh_atom = int(shell_atom) + 1
      sh_l = int(shell_l)
      sh_nprim = int(shell_nprim)
      if (any(sh_nprim < 1)) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_internal_api", &
                        "Every shell needs at least one primitive")
         return
      end if
      nprim_tot = sum(sh_nprim)
      call c_f_pointer(c_exps, exps, [nprim_tot])
      call c_f_pointer(c_coeffs, coeffs, [nprim_tot])

      p_scale => null()
      if (c_associated(c_scale)) then
         call c_f_pointer(c_scale, p_scale)
         use_scale = p_scale
      else
         use_scale = 1000.0_wp
      end if

      p_debug => null()
      if (c_associated(c_debug)) then
         call c_f_pointer(c_debug, p_debug)
         use_debug = p_debug
      else
         use_debug = .false.
      end if

      p_verbose => null()
      if (c_associated(c_verbose)) then
         call c_f_pointer(c_verbose, p_verbose)
         use_verbose = p_verbose
      else
         use_verbose = 0
      end if

      p_do_fine => null()
      if (c_associated(c_do_fine)) then
         call c_f_pointer(c_do_fine, p_do_fine)
         use_do_fine = p_do_fine
      else
         use_do_fine = .false.
      end if

      p_wleb_prune_level => null()
      if (c_associated(c_wleb_prune_level)) then
         call c_f_pointer(c_wleb_prune_level, p_wleb_prune_level)
         use_wleb_prune_level = p_wleb_prune_level
      else
         use_wleb_prune_level = 0
      end if

      call new_radii("cpcm", radius_model, cavity_error)
      if (allocated(cavity_error)) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_internal_api", &
                        cavity_error%message)
         return
      end if

      ! Parse optional master tolerance (default: compiled DROP default)
      if (.not. decode_drop_tolerance(c_tolerance, use_tolerance)) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_internal_api", &
                        "DROP tolerance must be positive")
         return
      end if

      allocate (cav)
      allocate (cavity_type_drop :: cav%ptr)
      call new_context(cav%ctx, verbosity=use_verbose, debug=use_debug)
      pnleb => null()
      if (c_associated(nleb)) call c_f_pointer(nleb, pnleb)

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         block
            type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf_template

            call lsf_template%new(sh_atom, sh_l, sh_nprim, exps, coeffs, &
                                  real(rho_iso, wp), scale=use_scale, error=cavity_error)
            if (allocated(cavity_error)) then
               call api_error(error%ptr, "new_drop_cavity_isodensity_internal_api", &
                              cavity_error%message)
               if (associated(cav%ptr)) deallocate (cav%ptr)
               deallocate (cav)
               return
            end if

            if (associated(pnleb)) then
               call new_cavity_drop(cavity, cav%ctx, nleb=pnleb, &
                                    tolerance=use_tolerance, &
                                    do_fine=use_do_fine, &
                                    wleb_prune_level=use_wleb_prune_level, &
                                    radius_model=radius_model, lsf_model=lsf_template, &
                                    error=cavity_error)
            else
               call new_cavity_drop(cavity, cav%ctx, &
                                    tolerance=use_tolerance, &
                                    do_fine=use_do_fine, wleb_prune_level=use_wleb_prune_level, &
                                    radius_model=radius_model, &
                                    lsf_model=lsf_template, error=cavity_error)
            end if
         end block
      end select

      if (allocated(cavity_error)) then
         call api_error(error%ptr, "new_drop_cavity_isodensity_internal_api", &
                        cavity_error%message)
         if (associated(cav%ptr)) deallocate (cav%ptr)
         deallocate (cav)
         return
      end if

      vcav = c_loc(cav)
   end function new_drop_cavity_isodensity_internal_api

!> Report the internal isodensity cartesian-component layout.
!>
!> Two-pass: call once with the array pointers NULL to read ncart/nshell, then
!> allocate and call again.  shell_cart_offset is 0-based (length nshell+1);
!> comp_lx/ly/lz are the per-component monomial powers (length ncart) that define
!> the ordering the host must match when building its density transform.
   subroutine get_isodensity_cart_layout_api(verror, vcav, c_ncart, c_nshell, &
         c_shell_off, c_comp_lx, c_comp_ly, c_comp_lz) &
         & bind(C, name=namespace//"get_isodensity_cart_layout")
      type(c_ptr), value :: verror
      type(c_ptr), value :: vcav
      type(c_ptr), value :: c_ncart
      type(c_ptr), value :: c_nshell
      type(c_ptr), value :: c_shell_off
      type(c_ptr), value :: c_comp_lx
      type(c_ptr), value :: c_comp_ly
      type(c_ptr), value :: c_comp_lz

      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav
      integer(c_int), pointer :: p_ncart, p_nshell
      integer(c_int), pointer :: shell_off(:), comp_lx(:), comp_ly(:), comp_lz(:)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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

!> Install the cartesian-monomial density matrix for the internal isodensity LSF.
!>
!> dcart is the ncart-by-ncart density matrix in moist's cartesian-monomial basis
!> (column-major), typically the host's density transformed by the cart<-basis
!> map derived from get_isodensity_cart_layout.  Call before each cavity build.
   subroutine set_isodensity_density_api(verror, vcav, ncart, c_dcart) &
         & bind(C, name=namespace//"set_isodensity_density")
      type(c_ptr), value :: verror
      type(c_ptr), value :: vcav
      integer(c_int), value :: ncart
      type(c_ptr), value :: c_dcart

      type(vp_error), pointer :: error
      type(vp_cavity), pointer :: cav
      real(c_double), pointer :: dcart(:, :)
      type(error_type), allocatable :: cavity_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "set_isodensity_density", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)
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

!> Return the master numerical tolerance configured on a DROP cavity.
!> @param[in]  verror Opaque API error handle
!> @param[in]  vcav Opaque cavity handle
!> @param[out] c_tolerance Configured master tolerance
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

!> Rebuild DROP cavity for a new geometry
   subroutine update_drop_cavity_api(verror, vcav, vmol, nleb) &
         & bind(C, name=namespace//"update_drop_cavity")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(c_ptr), value :: vmol
      type(vp_structure), pointer :: mol
      type(c_ptr), value :: nleb
      integer(c_int), pointer :: pnleb
      type(error_type), allocatable :: cavity_error
      integer :: nat

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "update_drop_cavity_api", "DROP cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "update_drop_cavity_api", "Cavity is not initialized")
         return
      end if
      if (.not. cav%owned) then
         call api_error(error%ptr, "update_drop_cavity_api", &
                        "Cannot update a borrowed cavity handle")
         return
      end if

      if (.not. c_associated(vmol)) then
         call api_error(error%ptr, "update_drop_cavity_api", "Molecular structure data is missing")
         return
      end if
      call c_f_pointer(vmol, mol)

      nat = mol%ptr%nat
      if (nat <= 0) then
         call api_error(error%ptr, "update_drop_cavity_api", "Invalid number of atoms")
         return
      end if

      pnleb => null()
      if (c_associated(nleb)) call c_f_pointer(nleb, pnleb)

      ! Use SELECT TYPE to access DROP-specific functionality
      if (associated(pnleb)) then
         select type (cavity => cav%ptr)
         type is (cavity_type_drop)
            if (.not. allocated(cavity%radius_model)) then
               call api_error(error%ptr, "update_drop_cavity_api", "Cavity radius model is not initialized")
               return
            end if
            if (.not. allocated(cavity%lsf_model)) then
               call api_error(error%ptr, "update_drop_cavity_api", "Cavity LSF model is not initialized")
               return
            end if
            block
               !> Detached copies of everything the rebuild has to carry over.
               !> The constructor takes the cavity itself as `intent(inout)` and
               !> resets `param` (and deallocates the two models) before reading
               !> its own dummy arguments, so passing cavity components directly
               !> would hand it storage it has already overwritten or released.
               class(radius_type), allocatable :: radius_keep
               class(moist_cavity_drop_lsf_type), allocatable :: lsf_keep
               !> Configured numerical knobs to survive the rebuild
               real(wp) :: keep_tolerance
               integer :: keep_proj_maxiter, keep_proj_level, keep_wleb_prune_level

               allocate (radius_keep, source=cavity%radius_model)
               allocate (lsf_keep, source=cavity%lsf_model)
               keep_tolerance = cavity%param%tolerance
               keep_proj_maxiter = cavity%param%proj_maxiter
               keep_proj_level = cavity%param%proj_level
               keep_wleb_prune_level = cavity%param%wleb_prune_level

               !> Re-running the constructor resets every parameter to its
               !> compiled default, so carry the configured knobs across the
               !> rebuild; only the Lebedev order is meant to change here.
               call new_cavity_drop(cavity, ctx=cav%ctx, nleb=pnleb, &
                                    tolerance=keep_tolerance, &
                                    proj_maxiter=keep_proj_maxiter, &
                                    proj_level=keep_proj_level, &
                                    wleb_prune_level=keep_wleb_prune_level, &
                                    radius_model=radius_keep, lsf_model=lsf_keep, error=cavity_error)
            end block
         end select
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "update_drop_cavity_api", cavity_error%message)
            return
         end if
      end if

      call cav%ptr%update(mol%ptr, error=cavity_error)
      if (allocated(cavity_error)) then
         call api_error(error%ptr, "update_drop_cavity_api", cavity_error%message)
         return
      end if

   end subroutine update_drop_cavity_api

!> Get all size information (more efficient than multiple calls)
   subroutine get_drop_sizes_api(verror, vcav, ngrid, nmax, nsph) &
         & bind(C, name=namespace//"get_drop_sizes")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      integer(c_int), intent(out) :: ngrid
      integer(c_int), intent(out) :: nmax
      integer(c_int), intent(out) :: nsph

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_drop_sizes_api", "DROP cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_drop_sizes_api", "Cavity is not initialized")
         return
      end if

      ! Use SELECT TYPE for DROP-specific fields
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         ngrid = cavity%ngrid
         nmax = cavity%nmax
         nsph = cavity%nsph
      class default
         call api_error(error%ptr, "get_drop_sizes_api", "Cavity is not DROP type")
      end select

   end subroutine get_drop_sizes_api

!> Get grid sizes (total grid points and raw grid size before filtering)
   subroutine get_drop_grid_size_api(verror, vcav, ngrid, nmax) &
         & bind(C, name=namespace//"get_drop_grid_size")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      integer(c_int), intent(out) :: ngrid
      integer(c_int), intent(out) :: nmax

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_drop_grid_size_api", "DROP cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_drop_grid_size_api", "Cavity is not initialized")
         return
      end if

      ! Use SELECT TYPE for DROP-specific fields
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         ngrid = cavity%ngrid
         nmax = cavity%nmax
      class default
         call api_error(error%ptr, "get_drop_grid_size_api", "Cavity is not DROP type")
      end select

   end subroutine get_drop_grid_size_api

!> Get all DROP cavity results in one call
!> This is more efficient than calling individual getters
!> The caller passes the capacities it allocated the arrays with; they have to
!> be at least the cavity's own ngrid/nsph, otherwise nothing is written.
   subroutine get_drop_results_api(verror, vcav, ngrid_cap, nsph_cap, &
         & area, volume, ngrid, nmax, nsph, &
         & xyz, normal0, wleb, a, r_iI0, &
         & f, rho, &
         & owner, converged, radii, asph) &
         & bind(C, name=namespace//"get_drop_results")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav

      ! Caller-allocated array capacities (must be >= the cavity sizes)
      integer(c_int), value :: ngrid_cap
      integer(c_int), value :: nsph_cap

      ! Scalar outputs
      real(c_double), intent(out) :: area
      real(c_double), intent(out) :: volume
      integer(c_int), intent(out) :: ngrid
      integer(c_int), intent(out) :: nmax
      integer(c_int), intent(out) :: nsph

      ! Array outputs, dimensioned with the caller's capacity; only the leading
      ! ngrid (resp. nsph) entries of each are written.
      ! Grid point data (capacity: ngrid_cap)
      real(c_double), intent(out) :: xyz(3, ngrid_cap)
      real(c_double), intent(out) :: normal0(3, ngrid_cap)
      real(c_double), intent(out) :: wleb(ngrid_cap)
      real(c_double), intent(out) :: a(ngrid_cap)
      real(c_double), intent(out) :: r_iI0(ngrid_cap)
      real(c_double), intent(out) :: f(ngrid_cap)
      real(c_double), intent(out) :: rho(ngrid_cap)
      integer(c_int), intent(out) :: owner(ngrid_cap)
      logical(c_bool), intent(out) :: converged(ngrid_cap)

      ! Per-sphere data (capacity: nsph_cap)
      real(c_double), intent(out) :: radii(nsph_cap)
      real(c_double), intent(out) :: asph(nsph_cap)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_drop_results_api", "DROP cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_drop_results_api", "Cavity is not initialized")
         return
      end if

      if (.not. allocated(cav%ptr%total_area)) then
         call api_error(error%ptr, "get_drop_results_api", "DROP cavity is not initialized")
         return
      end if

      ! Use SELECT TYPE to access DROP-specific fields
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         if (ngrid_cap < cavity%ngrid .or. nsph_cap < cavity%nsph) then
            call api_error(error%ptr, "get_drop_results_api", &
                           "Array capacity too small - use values from get_cavity_sizes")
            return
         end if

         ! Get scalar values
         area = cavity%total_area
         volume = cavity%total_volume
         ngrid = cavity%ngrid
         nmax = cavity%nmax
         nsph = cavity%nsph

         ! Get grid point arrays
         xyz(:, :cavity%ngrid) = cavity%xyz(:, :)
         normal0(:, :cavity%ngrid) = cavity%normal0(:, :)
         wleb(:cavity%ngrid) = cavity%wleb
         a(:cavity%ngrid) = cavity%a
         r_iI0(:cavity%ngrid) = cavity%r_iI0
         f(:cavity%ngrid) = cavity%f
         rho(:cavity%ngrid) = cavity%rho
         ! Convert from Fortran 1-based to C 0-based indexing
         owner(:cavity%ngrid) = cavity%owner - 1
         converged(:cavity%ngrid) = cavity%converged

         ! Get per-sphere arrays
         radii(:cavity%nsph) = cavity%radii
         asph(:cavity%nsph) = cavity%asph
      class default
         call api_error(error%ptr, "get_drop_results_api", "Cavity is not DROP type")
      end select

   end subroutine get_drop_results_api

!> Assemble A-matrix and compute xi values
!> This must be called before accessing xi values or using the A-matrix
!> `ngrid_cap` is the caller's allocated grid capacity and has to be at least
!> the cavity's own ngrid.
   subroutine assemble_amat_api(verror, vcav, ngrid_cap, amat0, xi) &
         & bind(C, name=namespace//"assemble_amat")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      integer(c_int), value :: ngrid_cap
      real(c_double), intent(out) :: amat0(ngrid_cap, ngrid_cap)
      real(c_double), intent(out) :: xi(ngrid_cap)
      real(wp), allocatable :: amat0_local(:, :)
      type(error_type), allocatable :: cavity_error
      integer :: ngrid

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "assemble_amat_api", "DROP cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "assemble_amat_api", "DROP cavity is not initialized")
         return
      end if

      associate (cavity => cav%ptr)
         ngrid = cavity%ngrid
         if (ngrid_cap < ngrid) then
            call api_error(error%ptr, "assemble_amat_api", &
                           "Array capacity too small - use ngrid value from get_cavity_sizes")
            return
         end if
         if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f) &
             .or. .not. allocated(cavity%xyz)) then
            call api_error(error%ptr, "assemble_amat_api", &
                           "Cavity does not expose Gaussian surface widths")
            return
         end if

         ! Generic Gaussian surface-charge interaction matrix
         allocate (amat0_local(ngrid, ngrid))
         call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, amat0_local, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "assemble_amat_api", cavity_error%message)
            return
         end if

         ! Copy results to output arrays
         amat0(:ngrid, :ngrid) = amat0_local(:, :)
         xi(:cavity%ngrid) = cavity%xi0
      end associate

   end subroutine assemble_amat_api

!> Return Gaussian widths and switching factors for a Gaussian PCM cavity.
   subroutine get_cavity_gaussian_api(verror, vcav, ngrid_cap, c_xi, c_f) &
         & bind(C, name=namespace//"get_cavity_gaussian")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      integer(c_int), value :: ngrid_cap
      type(c_ptr), value :: c_xi, c_f
      real(c_double), pointer :: xi(:), f(:)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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
                        "Array capacity too small - use ngrid value from get_cavity_sizes")
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

!===============================================================================!
! GENERIC CAVITY API (Tier 1 - works on all cavity types)
!===============================================================================!

!> Generic update cavity - works for all cavity types
   subroutine update_cavity_api(verror, vcav, vmol) &
         & bind(C, name=namespace//"update_cavity")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(c_ptr), value :: vmol
      type(vp_structure), pointer :: mol
      type(error_type), allocatable :: cavity_error
      integer :: nat

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "update_cavity_api", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "update_cavity_api", "Cavity is not initialized")
         return
      end if
      if (.not. cav%owned) then
         call api_error(error%ptr, "update_cavity_api", &
                        "Cannot update a borrowed cavity handle")
         return
      end if

      if (.not. c_associated(vmol)) then
         call api_error(error%ptr, "update_cavity_api", "Molecular structure data is missing")
         return
      end if
      call c_f_pointer(vmol, mol)

      nat = mol%ptr%nat
      if (nat <= 0) then
         call api_error(error%ptr, "update_cavity_api", "Invalid number of atoms")
         return
      end if

      ! Call deferred procedure - works for all cavity types
      call cav%ptr%update(mol%ptr, error=cavity_error)
      if (allocated(cavity_error)) then
         call api_error(error%ptr, "update_cavity_api", cavity_error%message)
         return
      end if

   end subroutine update_cavity_api

!> Get generic cavity sizes - works for all cavity types
   subroutine get_cavity_sizes_api(verror, vcav, ngrid, nsph) &
         & bind(C, name=namespace//"get_cavity_sizes")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      integer(c_int), intent(out) :: ngrid
      integer(c_int), intent(out) :: nsph

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_cavity_sizes_api", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_cavity_sizes_api", "Cavity is not initialized")
         return
      end if

      if (.not. allocated(cav%ptr%radii)) then
         call api_error(error%ptr, "get_cavity_sizes_api", "Cavity is not built yet")
         return
      end if

      ! These fields are in base cavity_type
      ngrid = cav%ptr%ngrid
      nsph = size(cav%ptr%radii)

   end subroutine get_cavity_sizes_api

!> Get generic cavity results - works for all cavity types
!> Returns only fields from base cavity_type
!> The caller passes the capacities it allocated the arrays with; they have to
!> be at least the cavity's own ngrid/nsph, otherwise nothing is written.
   subroutine get_cavity_results_api(verror, vcav, ngrid_cap, nsph_cap, &
         & area, volume, ngrid, nsph, &
         & xyz, a, owner, converged, radii, asph) &
         & bind(C, name=namespace//"get_cavity_results")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav

      ! Caller-allocated array capacities (must be >= the cavity sizes)
      integer(c_int), value :: ngrid_cap
      integer(c_int), value :: nsph_cap

      ! Scalar outputs
      real(c_double), intent(out) :: area
      real(c_double), intent(out) :: volume
      integer(c_int), intent(out) :: ngrid
      integer(c_int), intent(out) :: nsph

      ! Array outputs, dimensioned with the caller's capacity; only the leading
      ! ngrid (resp. nsph) entries of each are written.
      ! Grid point data (capacity: ngrid_cap)
      real(c_double), intent(out) :: xyz(3, ngrid_cap)
      real(c_double), intent(out) :: a(ngrid_cap)
      integer(c_int), intent(out) :: owner(ngrid_cap)
      logical(c_bool), intent(out) :: converged(ngrid_cap)

      ! Per-sphere data (capacity: nsph_cap)
      real(c_double), intent(out) :: radii(nsph_cap)
      real(c_double), intent(out) :: asph(nsph_cap)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_cavity_results_api", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_cavity_results_api", "Cavity is not initialized")
         return
      end if

      if (.not. allocated(cav%ptr%total_area)) then
         call api_error(error%ptr, "get_cavity_results_api", "Cavity is not built yet")
         return
      end if

      if (ngrid_cap < cav%ptr%ngrid .or. nsph_cap < size(cav%ptr%radii)) then
         call api_error(error%ptr, "get_cavity_results_api", &
                        "Array capacity too small - use values from get_cavity_sizes")
         return
      end if

      ! Get scalar values from base cavity_type
      area = cav%ptr%total_area
      volume = cav%ptr%total_volume
      ngrid = cav%ptr%ngrid
      nsph = size(cav%ptr%radii)

      ! Get grid point arrays from base cavity_type
      xyz(:, :ngrid) = cav%ptr%xyz(:, :)
      a(:ngrid) = cav%ptr%a
      ! Convert from Fortran 1-based to C 0-based indexing
      owner(:ngrid) = cav%ptr%owner - 1
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         converged(:ngrid) = cavity%converged
      class default
         converged(:ngrid) = .true.
      end select

      ! Get per-sphere arrays from base cavity_type
      radii(:nsph) = cav%ptr%radii
      asph(:nsph) = cav%ptr%asph

   end subroutine get_cavity_results_api

!===============================================================================!
! TYPE-SPECIFIC CAVITY API (Tier 2 - DROP-specific fields)
!===============================================================================!

!> Get DROP-specific cavity data
!> Only works for cavity_type_drop, returns error for other types
! TODO: this routine should also accept NULL pointers
   subroutine get_drop_specific_api(verror, vcav, ngrid_cap, &
         & nmax, normal0, wleb, r_iI0, &
         & f, rho) &
         & bind(C, name=namespace//"get_drop_specific")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav

      ! Caller-allocated grid capacity (must be >= the cavity's ngrid)
      integer(c_int), value :: ngrid_cap

      ! DROP-specific outputs, dimensioned with the caller's capacity; only the
      ! leading ngrid entries of each are written.
      integer(c_int), intent(out) :: nmax
      real(c_double), intent(out) :: normal0(3, ngrid_cap)
      real(c_double), intent(out) :: wleb(ngrid_cap)
      real(c_double), intent(out) :: r_iI0(ngrid_cap)
      real(c_double), intent(out) :: f(ngrid_cap)
      real(c_double), intent(out) :: rho(ngrid_cap)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_drop_specific_api", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_drop_specific_api", "Cavity is not initialized")
         return
      end if

      ! Use SELECT TYPE to access DROP-specific fields
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         if (ngrid_cap < cavity%ngrid) then
            call api_error(error%ptr, "get_drop_specific_api", &
                           "Array capacity too small - use ngrid value from get_cavity_sizes")
            return
         end if

         ! Access DROP-specific fields
         nmax = cavity%nmax
         normal0(:, :cavity%ngrid) = cavity%normal0(:, :)
         wleb(:cavity%ngrid) = cavity%wleb
         r_iI0(:cavity%ngrid) = cavity%r_iI0
         f(:cavity%ngrid) = cavity%f
         rho(:cavity%ngrid) = cavity%rho
      class default
         call api_error(error%ptr, "get_drop_specific_api", "Cavity is not DROP type - cannot get DROP-specific data")
      end select

   end subroutine get_drop_specific_api

!> Get the stable per- numbering of a DROP cavity.
!> The numbering is a unique id (anchor id plus branch offset) that identifies
!> the same physical surface point across cavity rebuilds, even when points are
!> culled.  It lets callers match grid points between successive cavities (e.g.
!> to interpolate/animate the surface over an SCF or optimization trajectory).
!> @param[in]  verror     Error handle
!> @param[in]  vcav       DROP cavity handle
!> @param[in]  ngrid_cap  Caller-allocated capacity of numbering (>= ngrid)
!> @param[out] numbering  Per-point unique id, native point order (ngrid)
   subroutine get_drop_numbering_api(verror, vcav, ngrid_cap, numbering) &
         & bind(C, name=namespace//"get_drop_numbering")
      !> Error handle
      type(c_ptr), value :: verror
      !> Fortran error pointer
      type(vp_error), pointer :: error
      !> DROP cavity handle
      type(c_ptr), value :: vcav
      !> Fortran cavity pointer
      type(vp_cavity), pointer :: cav
      !> Caller-allocated grid capacity (must be >= the cavity's ngrid)
      integer(c_int), value :: ngrid_cap
      !> Per-point unique numbering; only the leading ngrid entries are written
      integer(c_int), intent(out) :: numbering(ngrid_cap)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_drop_numbering_api", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "get_drop_numbering_api", "Cavity is not initialized")
         return
      end if

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         if (.not. allocated(cavity%numbering)) then
            call api_error(error%ptr, "get_drop_numbering_api", &
               & "Cavity numbering is not available - update the cavity first")
            return
         end if
         if (ngrid_cap < cavity%ngrid) then
            call api_error(error%ptr, "get_drop_numbering_api", &
               & "Array capacity too small - use ngrid value from get_cavity_sizes")
            return
         end if
         numbering(:cavity%ngrid) = cavity%numbering
      class default
         call api_error(error%ptr, "get_drop_numbering_api", &
            & "Cavity is not DROP type - cannot get DROP-specific data")
      end select

   end subroutine get_drop_numbering_api

!===============================================================================!
! GRADIENT API (Tier 3 - Cavity and A-matrix gradients)
!===============================================================================!

!> Compute cavity gradient w.r.t. nuclear coordinates
!> Must be called after update_cavity and before get_cavity_gradient
   subroutine compute_cavity_gradient_api(verror, vcav) &
         & bind(C, name=namespace//"compute_cavity_gradient")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(error_type), allocatable :: cavity_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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
!> Must be called after update_cavity and before the *_rA contractions.
!> Restricts each grid point's nuclear coupling to its owner atom's rigid anchor
!> motion (the field's nuclear derivatives are zero for callback LSFs).
   subroutine compute_anchor_gradient_api(verror, vcav) &
         & bind(C, name=namespace//"compute_anchor_gradient")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(error_type), allocatable :: cavity_error

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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

!> Get the anchor-channel nuclear derivatives produced by the anchor pass.
!> Must call compute_anchor_gradient (or compute_cavity_gradient) first.
!> Arrays, written in Fortran shape notation (column-major: the leftmost index
!> is contiguous).  A C caller passes flat buffers of the same total size and
!> indexes element (i1,...,in) at i1 + d1*(i2 + d2*(...)), all 0-based:
!>   xyz1_rA(3, 3, nsph, ngrid)  - d(r_i)_j / d(R_A)_alpha  (j, alpha, A, grid)
!>   xi1_rA(3, nsph, ngrid)      - d(xi_i)  / d(R_A)_alpha   (alpha, A, grid)
!>   a_i1_rA(3, nsph, ngrid)     - d(a_i)   / d(R_A)_alpha   (alpha, A, grid)
!>   v_i1_rA(3, nsph, ngrid)     - d(v_i)   / d(R_A)_alpha   (alpha, A, grid)
!>   A_tot1_rA(3, nsph)          - d(total area)   / d(R_A)_alpha  (alpha, A)
!>   V_tot1_rA(3, nsph)          - d(total volume) / d(R_A)_alpha  (alpha, A)
!> The per-point area/volume elements (a_i1_rA, v_i1_rA) are the un-summed
!> counterparts of A_tot1_rA/V_tot1_rA (summing over the grid recovers the
!> totals).  The  area carries a switching-function dependence
!> (a_i ~ f_i / xi_i^2), so a_i1_rA is NOT recoverable from xi1_rA alone; it is
!> what a generic geometric surface functional (e.g. GOSTSHYP) needs for its
!> area route.
   subroutine get_anchor_gradient_api(verror, vcav, nsph_cap, ngrid_cap, xyz1_rA, xi1_rA, &
         & a_i1_rA, v_i1_rA, A_tot1_rA, V_tot1_rA) &
         & bind(C, name=namespace//"get_anchor_gradient")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      !> Caller-allocated capacities (must be >= the cavity's nsph/ngrid)
      integer(c_int), value :: nsph_cap
      integer(c_int), value :: ngrid_cap
      !> Output arrays, dimensioned with the caller's capacities; only the
      !> leading nsph/ngrid entries of each dimension are written.
      real(c_double), intent(out) :: xyz1_rA(3, 3, nsph_cap, ngrid_cap)
      real(c_double), intent(out) :: xi1_rA(3, nsph_cap, ngrid_cap)
      real(c_double), intent(out) :: a_i1_rA(3, nsph_cap, ngrid_cap)
      real(c_double), intent(out) :: v_i1_rA(3, nsph_cap, ngrid_cap)
      real(c_double), intent(out) :: A_tot1_rA(3, nsph_cap)
      real(c_double), intent(out) :: V_tot1_rA(3, nsph_cap)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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
                           "Array capacity too small - use values from get_cavity_sizes")
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
!> Must call compute_cavity_gradient first
!> Arrays, written in Fortran shape notation (column-major: the leftmost index
!> is contiguous).  A C caller passes flat buffers of the same total size:
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
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      !> Caller-allocated capacities (must be >= the cavity's nsph/ngrid)
      integer(c_int), value :: nsph_cap
      integer(c_int), value :: ngrid_cap

      ! Output arrays, dimensioned with the caller's capacities; only the
      ! leading nsph/ngrid entries of each dimension are written.
      real(c_double), intent(out) :: A_tot1_rA(3, nsph_cap)
      real(c_double), intent(out) :: V_tot1_rA(3, nsph_cap)
      real(c_double), intent(out) :: asph1_rA(3, nsph_cap, nsph_cap)
      real(c_double), intent(out) :: vsph1_rA(3, nsph_cap, nsph_cap)
      real(c_double), intent(out) :: xyz1_rA(3, 3, nsph_cap, ngrid_cap)
      real(c_double), intent(out) :: r_iI1_rA(3, nsph_cap, ngrid_cap)
      real(c_double), intent(out) :: rho1_rA(3, nsph_cap, ngrid_cap)

      integer :: nsph, ngrid

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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
                           "Array capacity too small - use values from get_cavity_sizes")
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
!> Must call compute_cavity_gradient first for gradient computation
!> Arrays, written in Fortran shape notation (column-major: the leftmost index
!> is contiguous).  A C caller passes flat buffers of the same total size:
!>   Amat0(ngrid, ngrid)                  - A-matrix (symmetric)
!>   Amat1_rA(3, nsph, ngrid, ngrid)      - gradient of A-matrix
!>   xi(ngrid)                            - xi values
   subroutine get_amat_gradient_api(verror, vcav, nsph_cap, ngrid_cap, &
         & Amat0, Amat1_rA, xi) &
         & bind(C, name=namespace//"get_amat_gradient")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      !> Caller-allocated capacities (must be >= the cavity's nsph/ngrid)
      integer(c_int), value :: nsph_cap
      integer(c_int), value :: ngrid_cap
      type(error_type), allocatable :: cavity_error

      ! Output arrays, dimensioned with the caller's capacities; only the
      ! leading nsph/ngrid entries of each dimension are written.
      real(c_double), intent(out) :: Amat0(ngrid_cap, ngrid_cap)
      real(c_double), intent(out) :: Amat1_rA(3, nsph_cap, ngrid_cap, ngrid_cap)
      real(c_double), intent(out) :: xi(ngrid_cap)

      real(wp), allocatable :: Amat0_f(:, :)
      real(wp), allocatable :: Amat1_rA_f(:, :, :, :)
      integer :: ngrid, nsph

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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
                           "Array capacity too small - use values from get_cavity_sizes")
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
!> Computes grad_rA = sum_ij q1_i (dA_ij/dR_A) q2_j
!> Must call compute_cavity_gradient first
   subroutine contract_amat1_q1q2_rA_api(verror, vcav, c_q1, c_q2, c_grad_rA) &
         & bind(C, name=namespace//"contract_amat1_q1q2_rA")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(c_ptr), value :: c_q1
      type(c_ptr), value :: c_q2
      type(c_ptr), value :: c_grad_rA
      real(c_double), pointer :: q1(:)
      real(c_double), pointer :: q2(:)
      real(c_double), pointer :: grad_rA(:, :)
      type(error_type), allocatable :: cavity_error
      integer :: ngrid, nsph

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

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
!> Computes weights w_xi, w_f, and w_xyz satisfying:
!>   q1^T dA q2 = sum_i w_xi_i dxi_i + w_f_i df_i + w_xyz(:,i).dxyz_i
   subroutine contract_amat1_q1q2_surface_weights_api(verror, vcav, c_q1, c_q2, &
         & c_w_xi, c_w_f, c_w_xyz) &
         & bind(C, name=namespace//"contract_amat1_q1q2_surface_weights")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(c_ptr), value :: c_q1
      type(c_ptr), value :: c_q2
      type(c_ptr), value :: c_w_xi
      type(c_ptr), value :: c_w_f
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

!> Contract the original DROP surface channels to LSF adjoint weights.
!> This wrapper preserves the version-0.5 ABI; use the extended entry point for
!> the normal and principal-curvature channels.
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

      call contract_surface_lsf_weights_extended_api(verror, vcav, c_w_xi, c_w_f, &
                                                     c_w_xyz, c_w_lsf0, c_w_lsf1, c_w_lsf2, &
                                                     c_null_ptr, c_null_ptr, c_null_ptr)

   end subroutine contract_surface_lsf_weights_api

!> Contract DROP surface weights to per-grid LSF adjoint weights.
!> The projected-coordinate and xi chains are always contracted; the optional
!> outward-normal (c_w_n) and principal-curvature (c_w_k1, c_w_k2) channels are
!> folded in when the caller supplies them (NULL skips the channel).
   subroutine contract_surface_lsf_weights_extended_api(verror, vcav, c_w_xi, c_w_f, c_w_xyz, &
         & c_w_lsf0, c_w_lsf1, c_w_lsf2, c_w_n, c_w_k1, c_w_k2) &
         & bind(C, name=namespace//"contract_surface_lsf_weights_extended")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(c_ptr), value :: c_w_xi
      type(c_ptr), value :: c_w_f
      type(c_ptr), value :: c_w_xyz
      type(c_ptr), value :: c_w_lsf0
      type(c_ptr), value :: c_w_lsf1
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

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "contract_surface_lsf_weights", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "contract_surface_lsf_weights", "Cavity is not initialized")
         return
      end if

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         ngrid = cavity%ngrid

         if (.not. c_associated(c_w_xi) .or. .not. c_associated(c_w_f) &
             .or. .not. c_associated(c_w_xyz) .or. .not. c_associated(c_w_lsf0) &
             .or. .not. c_associated(c_w_lsf1) .or. .not. c_associated(c_w_lsf2)) then
            call api_error(error%ptr, "contract_surface_lsf_weights", "Null array pointer provided")
            return
         end if

         call c_f_pointer(c_w_xi, w_xi, [ngrid])
         call c_f_pointer(c_w_f, w_f, [ngrid])
         call c_f_pointer(c_w_xyz, w_xyz, [3, ngrid])
         call c_f_pointer(c_w_lsf0, w_lsf0, [ngrid])
         call c_f_pointer(c_w_lsf1, w_lsf1, [3, ngrid])
         call c_f_pointer(c_w_lsf2, w_lsf2, [3, 3, ngrid])

         ! Optional channels: a NULL pointer leaves the matching accumulator
         ! channel at zero, so the contraction skips that channel entirely.
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
            call api_error(error%ptr, "contract_surface_lsf_weights", cavity_error%message)
            return
         end if

         call cavity%contract_surface_lsf_weights(acc, w_lsf0, w_lsf1, w_lsf2, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_surface_lsf_weights", cavity_error%message)
            return
         end if
      class default
         call api_error(error%ptr, "contract_surface_lsf_weights", &
                        "Cavity is not DROP type - contraction API only supports DROP cavities")
      end select

   end subroutine contract_surface_lsf_weights_extended_api

!> Contract nuclear + electronic CPCM terms with surface and electric-field data.
!> Must call compute_cavity_gradient first
   subroutine contract_nuc_elec_qefield_rA_api(verror, vcav, c_surface_q, c_qefield, c_za, c_grad_rA) &
         & bind(C, name=namespace//"contract_nuc_elec_qefield_rA")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      type(c_ptr), value :: c_surface_q
      type(c_ptr), value :: c_qefield
      type(c_ptr), value :: c_za
      type(c_ptr), value :: c_grad_rA
      real(c_double), pointer :: surface_q(:)
      real(c_double), pointer :: qefield(:, :)
      real(c_double), pointer :: za(:)
      real(c_double), pointer :: grad_rA(:, :)
      type(error_type), allocatable :: cavity_error
      integer :: ngrid, nsph

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "contract_nuc_elec_qefield_rA", "Cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. associated(cav%ptr)) then
         call api_error(error%ptr, "contract_nuc_elec_qefield_rA", "Cavity is not initialized")
         return
      end if

      associate (cavity => cav%ptr)
         if (.not. allocated(cavity%xyz) .or. .not. allocated(cavity%sphxyz) &
             .or. .not. allocated(cavity%xyz1_rA)) then
            call api_error(error%ptr, "contract_nuc_elec_qefield_rA", &
                           "Cavity position derivatives are unavailable - "// &
                           "call compute_cavity_gradient first")
            return
         end if

         ngrid = cavity%ngrid
         nsph = cavity%nsph

         if (.not. c_associated(c_surface_q) .or. .not. c_associated(c_qefield) &
             .or. .not. c_associated(c_za) .or. .not. c_associated(c_grad_rA)) then
            call api_error(error%ptr, "contract_nuc_elec_qefield_rA", "Null array pointer provided")
            return
         end if

         call c_f_pointer(c_surface_q, surface_q, [ngrid])
         call c_f_pointer(c_qefield, qefield, [3, ngrid])
         call c_f_pointer(c_za, za, [nsph])
         call c_f_pointer(c_grad_rA, grad_rA, [3, nsph])

         call pcm_electrostatic_nuclear_gradient(cavity%xyz, cavity%sphxyz, &
                                                 cavity%xyz1_rA, surface_q, qefield, za, &
                                                 grad_rA, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_nuc_elec_qefield_rA", cavity_error%message)
            return
         end if
      end associate

   end subroutine contract_nuc_elec_qefield_rA_api

!> Generic delete cavity - works for all cavity types
   subroutine delete_cavity_api(vcav) &
         & bind(C, name=namespace//"delete_cavity")
      type(c_ptr), intent(inout) :: vcav
      type(vp_cavity), pointer :: cav

      if (c_associated(vcav)) then
         call c_f_pointer(vcav, cav)
         if (cav%owned .and. associated(cav%ptr)) deallocate (cav%ptr)
         nullify (cav%ptr)
         call cav%ctx%delete()
         deallocate (cav)
         vcav = c_null_ptr
      end if

   end subroutine delete_cavity_api

!> Delete DROP cavity handle (legacy - use delete_cavity instead)
   subroutine delete_drop_cavity_api(vcav) &
         & bind(C, name=namespace//"delete_drop_cavity")
      type(c_ptr), intent(inout) :: vcav
      type(vp_cavity), pointer :: cav

      if (c_associated(vcav)) then
         call c_f_pointer(vcav, cav)
         if (cav%owned .and. associated(cav%ptr)) deallocate (cav%ptr)
         nullify (cav%ptr)
         call cav%ctx%delete()
         deallocate (cav)
         vcav = c_null_ptr
      end if

   end subroutine delete_drop_cavity_api

   subroutine f_c_character(rhs, lhs, len)
      character(kind=c_char), intent(out) :: lhs(*)
      character(len=*), intent(in) :: rhs
      integer, intent(in) :: len
      integer :: length

      if (len <= 0) return

      length = min(len - 1, len_trim(rhs))

      if (length > 0) lhs(1:length) = transfer(rhs(1:length), lhs(1:length))
      lhs(length + 1:length + 1) = c_null_char

   end subroutine f_c_character

   subroutine c_f_character(rhs, lhs, max_len)
      character(kind=c_char), intent(in) :: rhs(*)
      character(len=:, kind=c_char), allocatable, intent(out) :: lhs
      integer, intent(in) :: max_len
      integer :: ii

      do ii = 1, min(max_len, huge(ii) - 1)
         if (rhs(ii) == c_null_char) exit
      end do
      allocate (character(len=ii - 1) :: lhs)
      lhs = transfer(rhs(1:ii - 1), lhs)

   end subroutine c_f_character

   subroutine c_f_character_ptr(rhs_ptr, lhs, max_len, truncated)
      type(c_ptr), value, intent(in) :: rhs_ptr
      character(len=:, kind=c_char), allocatable, intent(out) :: lhs
      integer, intent(in) :: max_len
      logical, intent(out), optional :: truncated
      character(kind=c_char), pointer :: rhs(:)
      integer :: ii, nchar, scan_len
      logical :: has_null

      if (present(truncated)) truncated = .false.

      if (.not. c_associated(rhs_ptr)) then
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

      allocate (character(len=nchar, kind=c_char) :: lhs)
      if (nchar > 0) lhs = transfer(rhs(1:nchar), lhs)

   end subroutine c_f_character_ptr

!> Cold fusion check
   subroutine verify_structure(error, mol)
      type(error_type), allocatable, intent(out) :: error
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
