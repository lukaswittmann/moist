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
!    - Implement moist_new_iswig_cavity() constructor
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
   use moist_type, only: cavity_type, solvation_model, wavefunction_type
   use moist_radii, only: radius_type, new_radii, custom_radius_type
   use moist_radii_custom, only: new_custom_radii_atoms, new_custom_radii_elements
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_isodensity_callback, only: &
      moist_cavity_drop_lsf_isodensity_callback_type
   ! use moist_model_gems, only : gems_model, new_gems_model
   use moist_data_solvents, only: solvation_system_parameters, new_solvation_system_parameters, get_solvent_id
   use moist_version, only: get_moist_version
   use moist_output_ascii, only: ascii_moist_header => moist_header, &
      & HEADER_FULL, HEADER_SHORT, HEADER_ASCII, &
      ! & ascii_gems_header => gems_header, &
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
      class(cavity_type), pointer :: ptr => null()
      logical :: owned = .true.
   end type vp_cavity

   type :: vp_radii
      class(radius_type), allocatable :: ptr
   end type vp_radii

   type :: vp_model
      class(solvation_model), allocatable :: ptr
   end type vp_model

   public :: vp_error, vp_structure, vp_cavity, vp_radii, vp_model
   public :: get_version_api, get_version_string_api
   public :: new_error_api, check_error_api, check_error_exit_api, get_error_api, delete_error_api
   public :: new_structure_api, delete_structure_api, update_structure_api
   public :: new_cpcm_radii_api, new_smd_radii_api, new_d3_radii_api, new_cosmo_radii_api, new_bondi_radii_api
   public :: new_custom_radii_api, set_custom_radii_atoms_api, set_custom_radii_elements_api
   public :: delete_radii_api
   ! Solvation model API
   ! public :: new_gems_solvation_model_api
   ! public :: get_solvation_model_cavity_api
   public :: update_solvation_model_api
   public :: get_solvation_model_energy_api
   public :: delete_solvation_model_api
   ! Type-specific constructors
   public :: new_drop_cavity_api
   public :: new_drop_cavity_with_radii_api
   public :: new_drop_cavity_isodensity_callback_api
   ! Generic cavity operations
   public :: update_cavity_api
   public :: get_cavity_sizes_api
   public :: get_cavity_results_api
   public :: delete_cavity_api
   ! Type-specific getters
   public :: get_drop_specific_api
   public :: assemble_amat_api
   ! Gradient API
   public :: compute_cavity_gradient_api
   public :: get_cavity_gradient_api
   public :: get_amat_gradient_api
   public :: contract_amat1_q1q2_rA_api
   public :: contract_amat1_q1q2_surface_weights_api
   public :: contract_surface_lsf_weights_api
   public :: contract_nuc_elec_qefield_rA_api
   public :: print_header_api, print_header_api_short, print_header_api_ascii, print_version_api
   public :: print_build_header_api
   ! public :: print_gems_header_api

contains

!> API error helper - creates consistent error messages with routine context
   subroutine api_error(error, routine, msg)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: routine
      character(len=*), intent(in) :: msg
      call fatal_error(error, "["//routine//"] "//msg)
   end subroutine api_error

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
      write (unit, '(a,1x,i0,a,i0,a,i0)') "moist version", major, ".", minor, ".", patch

   end subroutine print_version_api

!> Print moist build information (version, git commit, compiler, host) to a file unit
   subroutine print_build_header_api(unit) &
         & bind(C, name=namespace//"print_build_header")
      integer(c_int), value, intent(in) :: unit
      call ascii_moist_build_header(unit)
   end subroutine print_build_header_api

! !> Print GEMS header banner to a file unit (use 6 for stdout)
! subroutine print_gems_header_api(unit) &
!       & bind(C, name=namespace//"print_gems_header")
!    integer(c_int), value, intent(in) :: unit
!    call ascii_gems_header(unit)
! end subroutine print_gems_header_api

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

!> Check error and exit with message if error is set
   subroutine check_error_exit_api(verror, charptr, buffersize) &
         & bind(C, name=namespace//"check_error_exit")
      use iso_fortran_env, only: error_unit
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      character(kind=c_char), intent(in) :: charptr(*)
      integer(c_int), intent(in), optional :: buffersize
      character(len=:), allocatable :: context
      integer :: max_length

      if (.not. c_associated(verror)) return

      call c_f_pointer(verror, error)

      if (.not. allocated(error%ptr)) return

      ! Convert C string to Fortran string for context
      if (present(buffersize)) then
         max_length = buffersize
      else
         max_length = 256
      end if
      call c_f_character(charptr, context, max_length)

      ! Print error message to stderr
      write (error_unit, '(a)') "[moist Error] "//trim(context)//": "//trim(error%ptr%message)

      ! Exit with failure status
      error stop 1

   end subroutine check_error_exit_api

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
      allocate (custom_radius_type :: radii%ptr)
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
      type is (custom_radius_type)
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
      type is (custom_radius_type)
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

! !> Create a new GEMS solvation model handle.
! !> The constructor is specific to GEMS; the returned handle uses the generic
! !> solvation model base type for later update/energy calls.
! function new_gems_solvation_model_api(verror, solvent, c_debug, c_verbose, c_read_parameters, &
!       & parameter_file) result(vmodel) &
!       & bind(C, name=namespace//"new_gems_solvation_model")
!    type(c_ptr), value :: verror
!    type(vp_error), pointer :: error
!    type(c_ptr), value :: solvent
!    type(c_ptr), value :: c_debug
!    type(c_ptr), value :: c_verbose
!    type(c_ptr), value :: c_read_parameters
!    type(c_ptr), value :: parameter_file
!    type(c_ptr) :: vmodel
!    type(vp_model), pointer :: model
!    type(solvation_system_parameters) :: system
!    character(len=:, kind=c_char), allocatable :: solvent_c
!    character(len=:), allocatable :: solvent_name
!    character(len=:, kind=c_char), allocatable :: parameter_file_c
!    integer :: solvent_id
!    logical(c_bool), pointer :: p_debug
!    integer(c_int), pointer :: p_verbose
!    logical(c_bool), pointer :: p_read_parameters
!    logical :: debug
!    integer :: verbosity
!    logical :: read_parameters
!    logical :: solvent_truncated
!    logical :: parameter_file_truncated
!    type(error_type), allocatable :: solvent_error
!    type(error_type), allocatable :: model_error

!    vmodel = c_null_ptr

!    if (.not.c_associated(verror)) return
!    call c_f_pointer(verror, error)

!    if (.not.c_associated(solvent)) then
!       call api_error(error%ptr, "new_gems_solvation_model_api", "Solvent name pointer is missing")
!       return
!    end if

!    call c_f_character_ptr(solvent, solvent_c, api_max_cstr, solvent_truncated)
!    if (len(solvent_c) == 0) then
!       call api_error(error%ptr, "new_gems_solvation_model_api", "Solvent name is empty")
!       return
!    end if
!    if (solvent_truncated) then
!       call api_error(error%ptr, "new_gems_solvation_model_api", &
!          & "Solvent name is too long or not null-terminated")
!       return
!    end if

!    allocate(character(len=len(solvent_c)) :: solvent_name)
!    solvent_name = transfer(solvent_c, solvent_name)

!    p_debug => null()
!    if (c_associated(c_debug)) then
!       call c_f_pointer(c_debug, p_debug)
!       debug = p_debug
!    else
!       debug = .false.
!    end if

!    p_verbose => null()
!    if (c_associated(c_verbose)) then
!       call c_f_pointer(c_verbose, p_verbose)
!       verbosity = p_verbose
!    else
!       verbosity = 0
!    end if

!    p_read_parameters => null()
!    if (c_associated(c_read_parameters)) then
!       call c_f_pointer(c_read_parameters, p_read_parameters)
!       read_parameters = p_read_parameters
!    else
!       read_parameters = .false.
!    end if

!    if (read_parameters) then
!       if (.not.c_associated(parameter_file)) then
!          call api_error(error%ptr, "new_gems_solvation_model_api", "parameter_file pointer is missing")
!          return
!       end if
!       call c_f_character_ptr(parameter_file, parameter_file_c, api_max_cstr, parameter_file_truncated)
!       if (parameter_file_truncated) then
!          call api_error(error%ptr, "new_gems_solvation_model_api", &
!             & "parameter_file is too long or not null-terminated")
!          return
!       end if
!       if (len(parameter_file_c) == 0) then
!          call api_error(error%ptr, "new_gems_solvation_model_api", "parameter_file is empty")
!          return
!       end if
!    else
!       allocate(character(len=0, kind=c_char) :: parameter_file_c)
!    end if

!    call get_solvent_id(solvent_name, solvent_id, solvent_error)
!    if (allocated(solvent_error)) then
!       call api_error(error%ptr, "new_gems_solvation_model_api", solvent_error%message)
!       return
!    end if
!    call new_solvation_system_parameters(system, solvent_id, error=solvent_error)
!    if (allocated(solvent_error)) then
!       call api_error(error%ptr, "new_gems_solvation_model_api", solvent_error%message)
!       return
!    end if

!    allocate(model)
!    allocate(gems_model :: model%ptr)

!    select type (gems => model%ptr)
!    type is (gems_model)
!       call new_gems_model(model_error, gems, system, debug=debug, verbosity=verbosity, &
!          read_parameters=read_parameters, parameter_file=parameter_file_c)
!    end select

!    if (allocated(model_error)) then
!       call api_error(error%ptr, "new_gems_solvation_model_api", model_error%message)
!       if (allocated(model%ptr)) deallocate(model%ptr)
!       deallocate(model)
!       return
!    end if

!    vmodel = c_loc(model)

! end function new_gems_solvation_model_api

!> Delete solvation model handle.
   subroutine delete_solvation_model_api(vmodel) &
         & bind(C, name=namespace//"delete_solvation_model")
      type(c_ptr), intent(inout) :: vmodel
      type(vp_model), pointer :: model

      if (c_associated(vmodel)) then
         call c_f_pointer(vmodel, model)

         if (allocated(model%ptr)) deallocate (model%ptr)
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

      call model%ptr%update(mol%ptr, error=model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "update_solvation_model_api", model_error%message)
      end if

   end subroutine update_solvation_model_api

!> Get the total solvation energy from a solvation model.
!> GEMS currently does not require additional wavefunction data for the API.
   subroutine get_solvation_model_energy_api(verror, vmodel, energy) &
         & bind(C, name=namespace//"get_solvation_model_energy")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vmodel
      type(vp_model), pointer :: model
      real(c_double), intent(out) :: energy
      type(wavefunction_type) :: wfn
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

      call model%ptr%get_energy(wfn, energy, model_error)
      if (allocated(model_error)) then
         call api_error(error%ptr, "get_solvation_model_energy_api", model_error%message)
      end if

   end subroutine get_solvation_model_energy_api

! !> Get a borrowed cavity handle from a solvation model.
! !> The returned handle is NOT owned by the caller - do NOT call moist_update_cavity
! !> or moist_delete_cavity on it. It remains valid as long as the parent model exists.
! !> Call moist_update_solvation_model before extracting the cavity.
! function get_solvation_model_cavity_api(verror, vmodel) result(vcav) &
!       & bind(C, name=namespace//"get_solvation_model_cavity")
!    type(c_ptr), value :: verror
!    type(vp_error), pointer :: error
!    type(c_ptr), value :: vmodel
!    type(vp_model), pointer :: model
!    type(c_ptr) :: vcav
!    type(vp_cavity), pointer :: cav

!    vcav = c_null_ptr

!    if (.not.c_associated(verror)) return
!    call c_f_pointer(verror, error)

!    if (.not.c_associated(vmodel)) then
!       call api_error(error%ptr, "get_solvation_model_cavity", "Model handle is missing")
!       return
!    end if
!    call c_f_pointer(vmodel, model)

!    if (.not.allocated(model%ptr)) then
!       call api_error(error%ptr, "get_solvation_model_cavity", "Model is not initialized")
!       return
!    end if

!    select type (gems => model%ptr)
!    type is (gems_model)
!       if (.not.allocated(gems%cavity)) then
!          call api_error(error%ptr, "get_solvation_model_cavity", &
!             "Model cavity is not initialized - call moist_update_solvation_model first")
!          return
!       end if
!       allocate(cav)
!       cav%ptr => gems%cavity
!       cav%owned = .false.
!       vcav = c_loc(cav)
!    class default
!       call api_error(error%ptr, "get_solvation_model_cavity", &
!          "This solvation model type does not expose a cavity")
!    end select

! end function get_solvation_model_cavity_api

!> Internal helper to create DROP cavity handles.
   subroutine new_drop_cavity_common(verror, vradii, nleb, c_debug, c_verbose, &
                                     c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine, routine_name, vcav)
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
      character(len=*), intent(in) :: routine_name
      type(c_ptr), intent(out) :: vcav
      type(vp_error), pointer :: error
      type(vp_radii), pointer :: radii_handle
      integer(c_int), pointer :: pnleb
      logical(c_bool), pointer :: p_debug
      integer(c_int), pointer :: p_verbose
      real(c_double), pointer :: p_blendk
      real(c_double), pointer :: p_blend1b
      real(c_double), pointer :: p_blend2b
      real(c_double), pointer :: p_blend3b
      logical(c_bool), pointer :: p_do_fine
      type(vp_cavity), pointer :: cav
      logical :: use_debug
      integer :: use_verbose
      real(wp) :: use_blendk, use_blend1b, use_blend2b, use_blend3b
      logical :: use_do_fine
      logical :: use_explicit_radii
      type(error_type), allocatable :: cavity_error
      class(radius_type), allocatable :: radius_model

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

      ! Parse optional blendk parameter (default: 2.0)
      p_blendk => null()
      if (c_associated(c_blendk)) then
         call c_f_pointer(c_blendk, p_blendk)
         use_blendk = p_blendk
      else
         use_blendk = 2.0_wp
      end if

      ! Parse optional blend1b parameter (default: 1.0)
      p_blend1b => null()
      if (c_associated(c_blend1b)) then
         call c_f_pointer(c_blend1b, p_blend1b)
         use_blend1b = p_blend1b
      else
         use_blend1b = 1.0_wp
      end if

      ! Parse optional blend2b parameter (default: 1.0)
      p_blend2b => null()
      if (c_associated(c_blend2b)) then
         call c_f_pointer(c_blend2b, p_blend2b)
         use_blend2b = p_blend2b
      else
         use_blend2b = 1.0_wp
      end if

      ! Parse optional blend3b parameter (default: 1.0)
      p_blend3b => null()
      if (c_associated(c_blend3b)) then
         call c_f_pointer(c_blend3b, p_blend3b)
         use_blend3b = p_blend3b
      else
         use_blend3b = 1.0_wp
      end if

      ! Parse optional do_fine flag (default: false)
      p_do_fine => null()
      if (c_associated(c_do_fine)) then
         call c_f_pointer(c_do_fine, p_do_fine)
         use_do_fine = p_do_fine
      else
         use_do_fine = .false.
      end if

      use_explicit_radii = c_associated(vradii)
      radii_handle => null()
      if (use_explicit_radii) then
         call c_f_pointer(vradii, radii_handle)
         if (.not. allocated(radii_handle%ptr)) then
            call api_error(error%ptr, routine_name, "Radii handle is not initialized")
            return
         end if
      end if

      allocate (cav)
      ! Allocate specific cavity_type_drop (not abstract base type)
      allocate (cavity_type_drop :: cav%ptr)
      pnleb => null()
      if (c_associated(nleb)) call c_f_pointer(nleb, pnleb)

      ! Use SELECT TYPE to access the DROP-specific type
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         if (.not. use_explicit_radii) then
            call new_radii("cpcm", radius_model, cavity_error)
            if (allocated(cavity_error)) then
               call api_error(error%ptr, routine_name, cavity_error%message)
               if (associated(cav%ptr)) deallocate (cav%ptr)
               deallocate (cav)
               return
            end if
         end if

         block
            type(moist_cavity_drop_lsf_svdw_type) :: svdw_template

            !> The cavity sets the LSF screening threshold from its own
            !> tolerance, so we only forward the shape parameters here.
            call svdw_template%new(blend_k=use_blendk, blend_1b=use_blend1b, &
                                   blend_2b=use_blend2b, blend_3b=use_blend3b)

            if (associated(pnleb)) then
               if (use_explicit_radii) then
                  call new_cavity_drop(cavity, nleb=pnleb, debug=use_debug, verbose=use_verbose, &
                                       do_fine=use_do_fine, radius_model=radii_handle%ptr, &
                                       lsf_model=svdw_template, error=cavity_error)
               else
                  call new_cavity_drop(cavity, nleb=pnleb, debug=use_debug, verbose=use_verbose, &
                                       do_fine=use_do_fine, radius_model=radius_model, &
                                       lsf_model=svdw_template, error=cavity_error)
               end if
            else
               if (use_explicit_radii) then
                  call new_cavity_drop(cavity, debug=use_debug, verbose=use_verbose, &
                                       do_fine=use_do_fine, radius_model=radii_handle%ptr, &
                                       lsf_model=svdw_template, error=cavity_error)
               else
                  call new_cavity_drop(cavity, debug=use_debug, verbose=use_verbose, &
                                       do_fine=use_do_fine, radius_model=radius_model, &
                                       lsf_model=svdw_template, error=cavity_error)
               end if
            end if
         end block
      end select
      if (allocated(cavity_error)) then
         call api_error(error%ptr, routine_name, cavity_error%message)
         if (associated(cav%ptr)) deallocate (cav%ptr)
         deallocate (cav)
         return
      end if

      ! Constructor only initializes - user must call update_cavity to build
      vcav = c_loc(cav)
   end subroutine new_drop_cavity_common

!> Create new DROP cavity handle with default CPCM radii model.
   function new_drop_cavity_api(verror, nleb, c_debug, c_verbose, &
         c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine) result(vcav) &
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
      type(c_ptr) :: vcav

      call new_drop_cavity_common(verror, c_null_ptr, nleb, c_debug, c_verbose, &
                                  c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine, "new_drop_cavity_api", vcav)
   end function new_drop_cavity_api

!> Create new DROP cavity handle with an explicit radii model.
   function new_drop_cavity_with_radii_api(verror, vradii, nleb, c_debug, c_verbose, &
         c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine) result(vcav) &
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
      type(c_ptr) :: vcav

      call new_drop_cavity_common(verror, vradii, nleb, c_debug, c_verbose, &
                                  c_blendk, c_blend1b, c_blend2b, c_blend3b, c_do_fine, "new_drop_cavity_with_radii_api", vcav)
   end function new_drop_cavity_with_radii_api

!> Create new DROP cavity handle using an external isodensity LSF callback.
   function new_drop_cavity_isodensity_callback_api(verror, callback, context, &
         c_scale, nleb, c_debug, c_verbose, c_do_fine, c_wleb_prune_level) result(vcav) &
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

      allocate (cav)
      allocate (cavity_type_drop :: cav%ptr)
      pnleb => null()
      if (c_associated(nleb)) call c_f_pointer(nleb, pnleb)

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         block
            type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf_template

            call lsf_template%new(callback, context, scale=use_scale)

            if (associated(pnleb)) then
               call new_cavity_drop(cavity, nleb=pnleb, debug=use_debug, &
                                    verbose=use_verbose, do_fine=use_do_fine, &
                                    wleb_prune_level=use_wleb_prune_level, &
                                    radius_model=radius_model, lsf_model=lsf_template, &
                                    error=cavity_error)
            else
               call new_cavity_drop(cavity, debug=use_debug, verbose=use_verbose, &
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
            call new_cavity_drop(cavity, nleb=pnleb, debug=.false., &
                                 radius_model=cavity%radius_model, lsf_model=cavity%lsf_model, error=cavity_error)
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
   subroutine get_drop_results_api(verror, vcav, &
         & area, volume, ngrid, nmax, nsph, &
         & xyz, normal0, wleb, a, r_iI0, &
         & f, rho, &
         & owner, converged, radii, asph) &
         & bind(C, name=namespace//"get_drop_results")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav

      ! Scalar outputs
      real(c_double), intent(out) :: area
      real(c_double), intent(out) :: volume
      integer(c_int), intent(out) :: ngrid
      integer(c_int), intent(out) :: nmax
      integer(c_int), intent(out) :: nsph

      ! Array outputs (caller must allocate with correct sizes)
      ! Grid point data (size: ngrid)
      real(c_double), intent(out) :: xyz(3, *)
      real(c_double), intent(out) :: normal0(3, *)
      real(c_double), intent(out) :: wleb(*)
      real(c_double), intent(out) :: a(*)
      real(c_double), intent(out) :: r_iI0(*)
      real(c_double), intent(out) :: f(*)
      real(c_double), intent(out) :: rho(*)
      integer(c_int), intent(out) :: owner(*)
      logical(c_bool), intent(out) :: converged(*)

      ! Per-sphere data (size: nsph)
      real(c_double), intent(out) :: radii(*)
      real(c_double), intent(out) :: asph(*)

      if (.not. c_associated(verror)) return
      call c_f_pointer(verror, error)

      if (.not. c_associated(vcav)) then
         call api_error(error%ptr, "get_drop_results_api", "DROP cavity handle is missing")
         return
      end if
      call c_f_pointer(vcav, cav)

      if (.not. allocated(cav%ptr%total_area)) then
         call api_error(error%ptr, "get_drop_results_api", "DROP cavity is not initialized")
         return
      end if

      ! Use SELECT TYPE to access DROP-specific fields
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
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
   subroutine assemble_amat_api(verror, vcav, ngrid, amat0, xi) &
         & bind(C, name=namespace//"assemble_amat")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      integer(c_int), intent(in) :: ngrid
      real(c_double), intent(out) :: amat0(ngrid, *)
      real(c_double), intent(out) :: xi(*)
      real(wp), allocatable :: amat0_local(:, :)
      type(error_type), allocatable :: cavity_error

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

      ! Use SELECT TYPE for DROP-specific functionality
      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         if (ngrid /= cavity%ngrid) then
            call api_error(error%ptr, "assemble_amat_api", &
                           "Size mismatch - use ngrid value from get_cavity_sizes")
            return
         end if

         ! Call the cavity's assemble_Amat012_rA method
         call cavity%Amat012_rA(amat0_local, error=cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "assemble_amat_api", cavity_error%message)
            return
         end if

         ! Copy results to output arrays
         amat0(:ngrid, :ngrid) = amat0_local(:, :)
         xi(:cavity%ngrid) = cavity%xi0
      class default
         call api_error(error%ptr, "assemble_amat_api", "Cavity is not DROP type - cannot assemble A-matrix")
      end select

   end subroutine assemble_amat_api

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

      ! These fields are in base cavity_type
      ngrid = cav%ptr%ngrid
      nsph = size(cav%ptr%radii)

   end subroutine get_cavity_sizes_api

!> Get generic cavity results - works for all cavity types
!> Returns only fields from base cavity_type
   subroutine get_cavity_results_api(verror, vcav, &
         & area, volume, ngrid, nsph, &
         & xyz, a, owner, converged, radii, asph) &
         & bind(C, name=namespace//"get_cavity_results")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav

      ! Scalar outputs
      real(c_double), intent(out) :: area
      real(c_double), intent(out) :: volume
      integer(c_int), intent(out) :: ngrid
      integer(c_int), intent(out) :: nsph

      ! Array outputs (caller must allocate with correct sizes)
      ! Grid point data (size: ngrid)
      real(c_double), intent(out) :: xyz(3, *)
      real(c_double), intent(out) :: a(*)
      integer(c_int), intent(out) :: owner(*)
      logical(c_bool), intent(out) :: converged(*)

      ! Per-sphere data (size: nsph)
      real(c_double), intent(out) :: radii(*)
      real(c_double), intent(out) :: asph(*)

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
      converged(:ngrid) = cav%ptr%converged

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
   subroutine get_drop_specific_api(verror, vcav, &
         & nmax, normal0, wleb, r_iI0, &
         & f, rho) &
         & bind(C, name=namespace//"get_drop_specific")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav

      ! DROP-specific outputs
      integer(c_int), intent(out) :: nmax
      real(c_double), intent(out) :: normal0(3, *)
      real(c_double), intent(out) :: wleb(*)
      real(c_double), intent(out) :: r_iI0(*)
      real(c_double), intent(out) :: f(*)
      real(c_double), intent(out) :: rho(*)

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
      call cav%ptr%get_gradient()

   end subroutine compute_cavity_gradient_api

!> Get cavity gradient arrays (DROP-specific)
!> Must call compute_cavity_gradient first
!> Arrays:
!>   A_tot1_rA(3, nsph)           - gradient of total area
!>   V_tot1_rA(3, nsph)           - gradient of total volume
!>   asph1_rA(3, nsph, nsph)      - gradient of per-sphere areas
!>   vsph1_rA(3, nsph, nsph)      - gradient of per-sphere volumes
!>   xyz1_rA(3, 3, nsph, ngrid) - grid point position derivatives (j, alpha, A, grid)
!>   r_iI1_rA(3, nsph, ngrid)     - gradient of grid-owner distances
!>   rho1_rA(3, nsph, ngrid)      - gradient of rho values
   subroutine get_cavity_gradient_api(verror, vcav, nsph_in, ngrid_in, &
         & A_tot1_rA, V_tot1_rA, asph1_rA, vsph1_rA, &
         & xyz1_rA, r_iI1_rA, rho1_rA) &
         & bind(C, name=namespace//"get_cavity_gradient")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      integer(c_int), intent(in) :: nsph_in
      integer(c_int), intent(in) :: ngrid_in

      ! Output arrays (caller must allocate with correct sizes from get_cavity_sizes)
      real(c_double), intent(out) :: A_tot1_rA(3, nsph_in)
      real(c_double), intent(out) :: V_tot1_rA(3, nsph_in)
      real(c_double), intent(out) :: asph1_rA(3, nsph_in, nsph_in)
      real(c_double), intent(out) :: vsph1_rA(3, nsph_in, nsph_in)
      real(c_double), intent(out) :: xyz1_rA(3, 3, nsph_in, ngrid_in)
      real(c_double), intent(out) :: r_iI1_rA(3, nsph_in, ngrid_in)
      real(c_double), intent(out) :: rho1_rA(3, nsph_in, ngrid_in)

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

         ! Validate input sizes
         if (nsph_in /= nsph .or. ngrid_in /= ngrid) then
            call api_error(error%ptr, "get_cavity_gradient", &
                           "Size mismatch - use values from get_cavity_sizes")
            return
         end if

         ! Copy gradient arrays
         A_tot1_rA(:, :) = cavity%A_tot1_rA(:, :)
         V_tot1_rA(:, :) = cavity%V_tot1_rA(:, :)
         asph1_rA(:, :, :) = cavity%asph1_rA(:, :, :)
         vsph1_rA(:, :, :) = cavity%vsph1_rA(:, :, :)
         xyz1_rA(:, :, :, :) = cavity%xyz1_rA(:, :, :, :)
         r_iI1_rA(:, :, :) = cavity%r_iI1_rA(:, :, :)
         rho1_rA(:, :, :) = cavity%rho1_rA(:, :, :)

      class default
         call api_error(error%ptr, "get_cavity_gradient", &
                        "Cavity is not DROP type - gradient API only supports DROP cavities")
      end select

   end subroutine get_cavity_gradient_api

!> Assemble A-matrix with gradient (DROP-specific)
!> Must call compute_cavity_gradient first for gradient computation
!> Arrays:
!>   Amat0(ngrid, ngrid)                  - A-matrix
!>   Amat1_rA(3, nsph, ngrid, ngrid)      - gradient of A-matrix
!>   xi(ngrid)                            - xi values
   subroutine get_amat_gradient_api(verror, vcav, nsph_in, ngrid_in, &
         & Amat0, Amat1_rA, xi) &
         & bind(C, name=namespace//"get_amat_gradient")
      type(c_ptr), value :: verror
      type(vp_error), pointer :: error
      type(c_ptr), value :: vcav
      type(vp_cavity), pointer :: cav
      integer(c_int), intent(in) :: nsph_in
      integer(c_int), intent(in) :: ngrid_in
      type(error_type), allocatable :: cavity_error

      ! Output arrays
      real(c_double), intent(out) :: Amat0(ngrid_in, ngrid_in)
      real(c_double), intent(out) :: Amat1_rA(3, nsph_in, ngrid_in, ngrid_in)
      real(c_double), intent(out) :: xi(ngrid_in)

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

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         ! Check if gradient was computed (needed for r_iI1_rA, f1_rA etc.)
         if (.not. allocated(cavity%r_iI1_rA)) then
            call api_error(error%ptr, "get_amat_gradient", &
                           "Gradient not computed - call compute_cavity_gradient first")
            return
         end if

         ngrid = cavity%ngrid
         nsph = cavity%nsph

         if (ngrid_in /= ngrid .or. nsph_in /= nsph) then
            call api_error(error%ptr, "get_amat_gradient", &
                           "Size mismatch - use values from get_cavity_sizes")
            return
         end if

         ! Call Amat012_rA with gradient
         call cavity%Amat012_rA(Amat0_f, Amat1_rA_f, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "get_amat_gradient", cavity_error%message)
            return
         end if

         ! Copy results to C arrays
         Amat0(:, :) = Amat0_f(:, :)
         Amat1_rA(:, :, :, :) = Amat1_rA_f(:, :, :, :)
         xi(:) = cavity%xi0(:)

      class default
         call api_error(error%ptr, "get_amat_gradient", &
                        "Cavity is not DROP type - A-matrix gradient only supports DROP cavities")
      end select

   end subroutine get_amat_gradient_api

!> Contract first derivatives of A with two grid vectors (DROP-specific)
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

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
         ngrid = cavity%ngrid
         nsph = cavity%nsph

         if (.not. c_associated(c_q1) .or. .not. c_associated(c_q2) .or. .not. c_associated(c_grad_rA)) then
            call api_error(error%ptr, "contract_amat1_q1q2_rA", "Null array pointer provided")
            return
         end if

         call c_f_pointer(c_q1, q1, [ngrid])
         call c_f_pointer(c_q2, q2, [ngrid])
         call c_f_pointer(c_grad_rA, grad_rA, [3, nsph])

         call cavity%contract_amat1_q1q2_rA(q1, q2, grad_rA, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_amat1_q1q2_rA", cavity_error%message)
            return
         end if
      class default
         call api_error(error%ptr, "contract_amat1_q1q2_rA", &
                        "Cavity is not DROP type - contraction API only supports DROP cavities")
      end select

   end subroutine contract_amat1_q1q2_rA_api

!> Contract first derivatives of A to per-grid surface weights (DROP-specific)
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

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
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

         call cavity%contract_amat1_q1q2_surface_weights(q1, q2, w_xi, w_f, w_xyz, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_amat1_q1q2_surface_weights", cavity_error%message)
            return
         end if
      class default
         call api_error(error%ptr, "contract_amat1_q1q2_surface_weights", &
                        "Cavity is not DROP type - contraction API only supports DROP cavities")
      end select

   end subroutine contract_amat1_q1q2_surface_weights_api

!> Contract DROP surface weights to per-grid LSF adjoint weights.
!> The current implementation includes the projected-coordinate chain only.
   subroutine contract_surface_lsf_weights_api(verror, vcav, c_w_xi, c_w_f, c_w_xyz, &
         & c_w_lsf0, c_w_lsf1, c_w_lsf2) &
         & bind(C, name=namespace//"contract_surface_lsf_weights")
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
      real(c_double), pointer :: w_xi(:)
      real(c_double), pointer :: w_f(:)
      real(c_double), pointer :: w_xyz(:, :)
      real(c_double), pointer :: w_lsf0(:)
      real(c_double), pointer :: w_lsf1(:, :)
      real(c_double), pointer :: w_lsf2(:, :, :)
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

         call cavity%contract_surface_lsf_weights(w_xi, w_f, w_xyz, &
                                                  w_lsf0, w_lsf1, w_lsf2, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_surface_lsf_weights", cavity_error%message)
            return
         end if
      class default
         call api_error(error%ptr, "contract_surface_lsf_weights", &
                        "Cavity is not DROP type - contraction API only supports DROP cavities")
      end select

   end subroutine contract_surface_lsf_weights_api

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

      select type (cavity => cav%ptr)
      type is (cavity_type_drop)
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

         call cavity%contract_nuc_elec_qefield_rA(surface_q, qefield, za, grad_rA, cavity_error)
         if (allocated(cavity_error)) then
            call api_error(error%ptr, "contract_nuc_elec_qefield_rA", cavity_error%message)
            return
         end if
      class default
         call api_error(error%ptr, "contract_nuc_elec_qefield_rA", &
                        "Cavity is not DROP type - contraction API only supports DROP cavities")
      end select

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
