!> Assemble standalone RISM runs from application settings and native objects
module moist_setup
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_input, only: rism_input_type, rism_grid_parameters_type, rism_potential_parameters_type
   use moist_model_type, only: solvation_model_type
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent
   use moist_model_moz_solvent_water, only: new_spce_water_geometry, is_water
   use moist_data_solvents, only: solvation_system_type, new_solvation_system, get_solvent_id
   use moist_model_moz_potential, only: moz_potential_type, lj_uff, lj_dreiding, lj_tm, lj_gaff, lj_oplsaa, &
      & lj_spce, lj_cspce, coulomb_spce, coulomb_cspce, custom_lj_type, new_custom_lj, &
      & fixed_charges_type, new_fixed_charges, &
      & lj_mixing_lorentz_berthelot, lj_mixing_geometric
   use moist_model_moz_potential_sites, only: potential_sites_type
   use moist_model_moz_potential_electrostatic_solvent, only: solvent_charges_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, default_molecular_recipe, &
      & moist_math_grid_atomic_shell_constant_type, new_constant_shell_policy
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_handymod_type, new_handymod_mapping
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   implicit none(type, external)
   private
   public :: new_rism_model, evaluate_rism_model, numerical_rism_gradient

contains

   !> Build updated bulk and UV objects; the UV model owns the internal bulk solve
   !>
   !> @param[inout] input Validated run settings
   !> @param[in] ctx Borrowed context, outlives the returned model
   !> @param[in] mol Solute, used for atom-count and charge checks
   !> @param[out] model Updated model, with solute terms attached
   !> @param[out] error Invalid input, solvent or native construction failure
   subroutine new_rism_model(input, ctx, mol, model, error)
      type(rism_input_type), intent(inout) :: input
      type(moist_context_type), intent(in), target :: ctx
      type(structure_type), intent(in) :: mol
      class(solvation_model_type), allocatable, intent(out) :: model
      type(error_type), allocatable, intent(out) :: error
      type(solvation_system_type) :: system
      type(solvent_vv_type) :: solvent
      type(structure_type) :: solvent_mol
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_3d_type), allocatable :: grid
      type(model_moz_1d_type), allocatable :: uv1d
      type(model_moz_3d_type), allocatable :: uv3d
      integer :: solvent_id

      call input%validate(error)
      if (allocated(error)) return
      call get_solvent_id(trim(input%system%solvent), solvent_id, error)
      if (allocated(error)) return
      call new_solvation_system(system, solvent_id, temperature=input%system%temperature, &
         & pressure_si=input%system%pressure, error=error)
      if (allocated(error)) return
      if (input%solvent_potential%lj == "spce" .or. input%solvent_potential%lj == "cspce") then
         if (is_water(system%solv_mol, solvent_id)) call new_spce_water_geometry(system%solv_mol)
      end if
      call new_vv_solvent(solvent, ctx, system, error)
      if (allocated(error)) return
      solvent_mol = solvent%structure()
      call configure_potential(solvent%potential, input%solvent_potential, solvent_mol%nat, error)
      if (.not. allocated(error)) call solvent%update(error)
      if (.not. allocated(error)) call check_sites(solvent%potential, error)
      if (allocated(error)) return

      select type (params => input%model_parameters)
      type is (model_moz_1d_parameters_type)
         params%solver%verbosity = max(0, ctx%verbosity())
         params%solver%nthreads = ctx%get_num_threads()
         params%solvent%solver%verbosity = params%solver%verbosity
         params%solvent%solver%nthreads = params%solver%nthreads
         call new_uniform_radial_pair(rgrid, kgrid, input%grid%npts, input%grid%dr, error)
         if (allocated(error)) return
         allocate(uv1d)
         call new_moz_1d_model(uv1d, ctx, rgrid, kgrid, solvent, error, params)
         if (.not. allocated(error)) call configure_potential(uv1d%potential, input%potential, mol%nat, error, &
            & mixing=solvent%potential%mixing())
         if (.not. allocated(error)) call uv1d%update(mol, error)
         if (.not. allocated(error)) call check_sites(uv1d%potential, error, mol%charge)
         if (allocated(error)) return
         call move_alloc(uv1d, model)
      type is (model_moz_3d_parameters_type)
         params%solver%verbosity = max(0, ctx%verbosity())
         params%solver%nthreads = ctx%get_num_threads()
         params%solvent%solver%verbosity = params%solver%verbosity
         params%solvent%solver%nthreads = params%solver%nthreads
         call new_spatial_grid(input%grid, grid, error)
         if (allocated(error)) return
         allocate(uv3d)
         call new_moz_3d_model(uv3d, ctx, grid, solvent, error, params)
         if (.not. allocated(error)) call configure_potential(uv3d%potential, input%potential, mol%nat, error, &
            & mixing=solvent%potential%mixing())
         if (.not. allocated(error)) call uv3d%update(mol, error)
         if (.not. allocated(error)) call check_sites(uv3d%potential, error, mol%charge)
         if (allocated(error)) return
         call move_alloc(uv3d, model)
      class default
         call fatal_error(error, "Unsupported native RISM parameter implementation")
      end select
   end subroutine new_rism_model

   !> Require the requested LJ/charge suppliers to have actually covered this side
   !>
   !> @param[in] potential Updated, parameter-only potential
   !> @param[out] error Missing supplier or inconsistent solute charge
   !> @param[in] total_charge Optional solute total charge, e
   subroutine check_sites(potential, error, total_charge)
      type(moz_potential_type), intent(in) :: potential
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: total_charge
      type(potential_sites_type) :: sites

      call potential%sites(sites, error)
      if (allocated(error)) return
      if (.not. allocated(sites%pair) .or. .not. allocated(sites%q)) then
         call fatal_error(error, "Selected LJ/charge models have no parameters for this structure")
         return
      end if
      if (present(total_charge)) then
         if (abs(sum(sites%q) - total_charge) > 1.0e-6_wp) then
            call fatal_error(error, "Solute charges must sum to the molecular charge (system.charge or --charge)")
         end if
      end if
   end subroutine check_sites

   !> Construct a Cartesian or bounded molecular point-grid template
   !>
   !> @param[in] params Spatial grid settings
   !> @param[out] grid Geometry-free grid, updated by the model
   !> @param[out] error Native grid or atomic-recipe construction failure
   subroutine new_spatial_grid(params, grid, error)
      type(rism_grid_parameters_type), intent(in) :: params
      class(moist_math_grid_3d_type), allocatable, intent(out) :: grid
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), allocatable :: cartesian
      type(moist_math_grid_3d_molecular_type), allocatable :: molecular
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_radial_mapping_handymod_type) :: mapping
      type(moist_math_grid_atomic_shell_constant_type) :: shells

      select case (params%kind)
      case ("cartesian")
         allocate(cartesian)
         call new_cartesian_point_grid(cartesian, error, nx=params%nx, ny=params%ny, nz=params%nz, &
            & dr=params%dr, margin=params%margin)
         if (.not. allocated(error)) call move_alloc(cartesian, grid)
      case ("molecular")
         call default_molecular_recipe(recipe, error)
         if (allocated(error)) return
         call new_handymod_mapping(mapping, 0.0_wp, params%rmax, params%handymod_exponent, error)
         if (.not. allocated(error)) call new_constant_shell_policy(shells, params%angular_degree, error)
         if (allocated(error)) return
         recipe%radial%npts = params%nrad
         deallocate(recipe%radial%mapping, recipe%shells)
         allocate(recipe%radial%mapping, source=mapping)
         allocate(recipe%shells, source=shells)
         allocate(molecular)
         call new_molecular_point_grid(molecular, error, recipe=recipe, dr=params%dr, &
            & kbuffer=params%kbuffer, nufft_tol=params%nufft_tolerance, becke_k=params%becke_k)
         if (.not. allocated(error)) call move_alloc(molecular, grid)
      case default
         call fatal_error(error, "Unknown spatial RISM grid: "//trim(params%kind))
      end select
   end subroutine new_spatial_grid

   !> Assemble parameter-only potential terms; the native owner updates them
   !>
   !> @param[inout] potential Empty potential
   !> @param[in] params Potential settings
   !> @param[in] natom Atom count of this side
   !> @param[out] error Array mismatch or native term/settings failure
   !> @param[in] mixing Optional common solvent mixing rule for UV tables
   subroutine configure_potential(potential, params, natom, error, mixing)
      type(moz_potential_type), intent(inout) :: potential
      type(rism_potential_parameters_type), intent(in) :: params
      integer, intent(in) :: natom
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: mixing
      type(fixed_charges_type) :: charges
      type(solvent_charges_type) :: table_charges
      real(wp), allocatable :: zero(:)
      integer :: rule

      call potential%set_ng_split(params%ng_split, params%alpha, error)
      if (allocated(error)) return
      rule = lj_mixing_lorentz_berthelot
      if (params%mixing == "geometric") rule = lj_mixing_geometric
      if (present(mixing)) rule = mixing
      call potential%set_mixing(rule, error)
      if (.not. allocated(error)) call add_lj(potential, params%lj, params, natom, error)
      if (.not. allocated(error) .and. len_trim(params%lj_fallback) > 0) then
         call add_lj(potential, params%lj_fallback, params, natom, error)
      end if
      if (allocated(error)) return
      select case (params%charges)
      case ("fixed")
         if (size(params%q) /= natom) then
            call fatal_error(error, "Potential q must contain one charge per atom, in structure order")
            return
         end if
         call new_fixed_charges(charges, params%q, error)
         if (.not. allocated(error)) call potential%add(charges, error)
      case ("none")
         allocate(zero(natom), source=0.0_wp)
         call new_fixed_charges(charges, zero, error)
         if (.not. allocated(error)) call potential%add(charges, error)
      case ("spce")
         call potential%add(coulomb_spce, error)
      case ("cspce")
         call potential%add(coulomb_cspce, error)
      case ("hirshfeld", "resp", "mbis", "chelpg")
         table_charges%model = params%charges
         table_charges%environment = params%environment
         call potential%add(table_charges, error)
      case default
         call fatal_error(error, "Unknown charge factory: "//trim(params%charges))
      end select
   end subroutine configure_potential

   !> Add one native LJ implementation to the potential's owned term list
   !>
   !> @param[inout] potential Potential being assembled
   !> @param[in] name Native LJ choice
   !> @param[in] params Explicit arrays for custom LJ
   !> @param[in] natom Atom count
   !> @param[out] error Unknown choice, array mismatch or native constructor failure
   subroutine add_lj(potential, name, params, natom, error)
      type(moz_potential_type), intent(inout) :: potential
      character(len=*), intent(in) :: name
      type(rism_potential_parameters_type), intent(in) :: params
      integer, intent(in) :: natom
      type(error_type), allocatable, intent(out) :: error
      type(custom_lj_type) :: custom

      select case (name)
      case ("uff")
         call potential%add(lj_uff, error)
      case ("dreiding")
         call potential%add(lj_dreiding, error)
      case ("tm")
         call potential%add(lj_tm, error)
      case ("gaff")
         call potential%add(lj_gaff, error)
      case ("oplsaa")
         call potential%add(lj_oplsaa, error)
      case ("spce")
         call potential%add(lj_spce, error)
      case ("cspce")
         call potential%add(lj_cspce, error)
      case ("custom")
         if (size(params%sigma) /= natom .or. size(params%epsilon) /= natom) then
            call fatal_error(error, "Custom LJ sigma and epsilon must contain one value per atom")
            return
         end if
         call new_custom_lj(custom, params%sigma, params%epsilon, error)
         if (.not. allocated(error)) call potential%add(custom, error)
      case default
         call fatal_error(error, "Unknown LJ factory: "//trim(name))
      end select
   end subroutine add_lj

   !> Evaluate a parameter-only model through the shared coupling lifecycle
   !>
   !> @param[inout] model Updated solvation model
   !> @param[out] energy Solvation free energy, Hartree
   !> @param[out] error Coupling, energy or gradient failure
   !> @param[out] gradient Optional nuclear gradient, Hartree/bohr
   subroutine evaluate_rism_model(model, energy, error, gradient)
      class(solvation_model_type), intent(inout), target :: model
      real(wp), intent(out) :: energy
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(out), optional :: gradient(:, :)
      type(coupling_type), pointer :: coupling
      type(response_type) :: response

      energy = 0.0_wp
      call model%new_coupling(coupling, error)
      if (allocated(error)) return
      call model%prepare_energy(coupling, error)
      if (.not. allocated(error)) call model%get_energy(coupling, energy, error)
      if (.not. allocated(error) .and. present(gradient)) then
         gradient = 0.0_wp
         call model%prepare_gradient(coupling, error)
         if (.not. allocated(error)) call model%get_gradient(coupling, response, gradient, error)
      end if
      call model%release_coupling(coupling)
   end subroutine evaluate_rism_model

   !> Central-difference CLI gradient with fixed charges, bulk and grid settings
   !>
   !> @param[inout] model Constructed model, restored to mol before returning
   !> @param[in] mol Original solute
   !> @param[out] gradient Nuclear gradient, Hartree/bohr
   !> @param[out] error Failed displaced or restored solve
   subroutine numerical_rism_gradient(model, mol, gradient, error)
      class(solvation_model_type), intent(inout), target :: model
      type(structure_type), intent(in) :: mol
      real(wp), intent(out) :: gradient(:, :)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: trial
      type(error_type), allocatable :: restore_error
      real(wp), parameter :: step = 1.0e-4_wp
      real(wp) :: plus, minus, restored
      integer :: iat, axis

      if (any(shape(gradient) /= [3, mol%nat])) then
         call fatal_error(error, "Numerical RISM gradient shape mismatch")
         return
      end if
      gradient = 0.0_wp
      atoms: do iat = 1, mol%nat
         do axis = 1, 3
            trial = mol
            trial%xyz(axis, iat) = mol%xyz(axis, iat) + step
            call model%update(trial, error)
            if (.not. allocated(error)) call evaluate_rism_model(model, plus, error)
            if (allocated(error)) exit atoms
            trial%xyz(axis, iat) = mol%xyz(axis, iat) - step
            call model%update(trial, error)
            if (.not. allocated(error)) call evaluate_rism_model(model, minus, error)
            if (allocated(error)) exit atoms
            gradient(axis, iat) = (plus - minus)/(2.0_wp*step)
         end do
      end do atoms
      call model%update(mol, restore_error)
      if (.not. allocated(error) .and. .not. allocated(restore_error)) then
         call evaluate_rism_model(model, restored, restore_error)
      end if
      if (.not. allocated(error)) call move_alloc(restore_error, error)
   end subroutine numerical_rism_gradient
end module moist_setup
