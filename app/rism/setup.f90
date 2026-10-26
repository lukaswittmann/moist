!> Assemble standalone RISM runs from application settings and native objects
module moist_app_rism_setup
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use mctc_io_utils, only: to_lower
   use moist_app_rism_input, only: rism_input_type, rism_checkpoint_files_type
   use moist_utils_env, only: get_env
   use moist_model_type, only: solvation_model_type
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent, solvent_vv_parameters_type
   use moist_model_moz_solvent_theory, only: canonical_bulk_theory
   use moist_model_moz_checkpoint_file, only: moz_file_checkpoint_type, new_file_checkpoint
   use moist_model_moz_solvent_water, only: new_spce_water_geometry, is_water
   use moist_data_solvents, only: get_solvent_id
   use moist_solvent_system, only: solvent_system_type, new_solvent_system
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_sites, only: potential_sites_type
   use moist_model_moz_potential_factory, only: configure_potential
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_factory, only: new_spatial_grid
   implicit none(type, external)
   private
   public :: new_rism_model, vv_checkpoint_name

   !> Closures whose automatic VV files can warm-start another closure, in order of preference
   character(len=4), parameter :: fallback_closures(5) = [character(len=4) :: "HNC", "PSE4", "PSE3", "PSE2", "KH"]

   !> Environment variable of the directory of automatic VV files
   character(len=*), parameter :: vv_outdir_env = "MOIST_RISM_VV_OUTDIR"

contains

   !> Build updated bulk and UV objects; the UV model owns the internal bulk solve
   !>
   !> - With `checkpoint.enabled`, the solvent and the UV model get file stores;
   !>   see attach_vv_checkpoint and uv_checkpoint
   !>
   !> @param[inout] input  validated run settings
   !> @param[in] ctx       borrowed context, outlives the returned model
   !> @param[in] mol       solute, used for atom-count and charge checks
   !> @param[out] model    updated model, with solute terms attached
   !> @param[out] error    invalid input, solvent or native construction failure
   subroutine new_rism_model(input, ctx, mol, model, error)
      !> Validated run settings
      type(rism_input_type), intent(inout) :: input
      !> Borrowed context, outlives the returned model
      type(moist_context_type), intent(in), target :: ctx
      !> Solute, used for atom-count and charge checks
      type(structure_type), intent(in) :: mol
      !> Updated model, with solute terms attached
      class(solvation_model_type), allocatable, intent(out) :: model
      !> Invalid input, solvent or native construction failure
      type(error_type), allocatable, intent(out) :: error
      type(solvent_system_type) :: system
      type(solvent_vv_type) :: solvent
      type(structure_type) :: solvent_mol
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_3d_type), allocatable :: grid
      type(model_moz_1d_type), allocatable :: uv1d
      type(model_moz_3d_type), allocatable :: uv3d
      type(moz_file_checkpoint_type) :: uv_store
      integer :: solvent_id

      call input%validate(error)
      if (allocated(error)) return
      call get_solvent_id(trim(input%system%solvent), solvent_id, error)
      if (allocated(error)) return
      call new_solvent_system(system, solvent_id, temperature=input%system%temperature, &
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
         if (input%checkpoint%enabled) then
            call attach_vv_checkpoint(solvent, input%checkpoint%vv, solvent_id, params%solvent)
            call uv_checkpoint(input%checkpoint%uv, "uv_rism_1d.mbin", uv_store)
         end if
         call new_uniform_radial_pair(rgrid, kgrid, input%grid%npts, input%grid%dr, error)
         if (allocated(error)) return
         allocate(uv1d)
         call new_moz_1d_model(uv1d, ctx, rgrid, kgrid, solvent, error, params)
         if (.not. allocated(error)) call configure_potential(uv1d%potential, input%potential, mol%nat, error, &
            & mixing=solvent%potential%mixing())
         if (.not. allocated(error)) call uv1d%update(mol, error)
         if (.not. allocated(error)) call check_sites(uv1d%potential, error, mol%charge)
         if (allocated(error)) return
         if (input%checkpoint%enabled) call uv1d%set_checkpoint(uv_store)
         call move_alloc(uv1d, model)
      type is (model_moz_3d_parameters_type)
         if (input%checkpoint%enabled) then
            call attach_vv_checkpoint(solvent, input%checkpoint%vv, solvent_id, params%solvent)
            call uv_checkpoint(input%checkpoint%uv, "uv_rism_3d.mbin", uv_store)
         end if
         call new_spatial_grid(input%grid, grid, error)
         if (allocated(error)) return
         allocate(uv3d)
         call new_moz_3d_model(uv3d, ctx, grid, solvent, error, params)
         if (.not. allocated(error)) call configure_potential(uv3d%potential, input%potential, mol%nat, error, &
            & mixing=solvent%potential%mixing())
         if (.not. allocated(error)) call uv3d%update(mol, error)
         if (.not. allocated(error)) call check_sites(uv3d%potential, error, mol%charge)
         if (allocated(error)) return
         if (input%checkpoint%enabled) call uv3d%set_checkpoint(uv_store)
         call move_alloc(uv3d, model)
      class default
         call fatal_error(error, "Unsupported native RISM parameter implementation")
      end select
   end subroutine new_rism_model

   !> Attach a VV file store to the solvent before the UV model copies it
   !>
   !> - Saved file: `checkpoint.vv.file`, else the automatic name (vv_checkpoint_name)
   !> - Guess files: `guess`, else `file`; with neither, the automatic files of the
   !>   requested closure and then of the other closures, which can only warm-start
   !>
   !> @param[inout] solvent  updated solvent
   !> @param[in]    files    VV checkpoint settings
   !> @param[in]    solvent_id  solvent id; 0 for a custom solvent
   !> @param[in]    params   settings of the VV solve
   subroutine attach_vv_checkpoint(solvent, files, solvent_id, params)
      !> Updated solvent
      type(solvent_vv_type), intent(inout) :: solvent
      !> VV checkpoint settings
      type(rism_checkpoint_files_type), intent(in) :: files
      !> Solvent id
      integer, intent(in) :: solvent_id
      !> Settings of the VV solve
      type(solvent_vv_parameters_type), intent(in) :: params
      type(moz_file_checkpoint_type) :: store
      character(len=:), allocatable :: file, guesses(:)
      character(len=:), allocatable :: theory
      integer :: i, n, width

      theory = canonical_bulk_theory(params%theory)
      if (allocated(files%file)) then
         file = trim(files%file)
      else
         file = vv_checkpoint_name(solvent_id, theory, params%closure)
      end if
      if (allocated(files%guess)) then
         guesses = [character(len=len_trim(files%guess)) :: trim(files%guess)]
      else if (allocated(files%file)) then
         guesses = [character(len=len(file)) :: file]
      else
         width = len(vv_checkpoint_name(solvent_id, theory, "PSE4")) + len_trim(params%closure)
         allocate (character(len=width) :: guesses(1 + size(fallback_closures)))
         guesses(1) = file
         n = 1
         do i = 1, size(fallback_closures)
            if (to_lower(trim(fallback_closures(i))) == to_lower(trim(adjustl(params%closure)))) cycle
            n = n + 1
            guesses(n) = vv_checkpoint_name(solvent_id, theory, fallback_closures(i))
         end do
         guesses = guesses(:n)
      end if
      call new_file_checkpoint(store, file, guesses, files%interval)
      call solvent%set_checkpoint(store)
   end subroutine attach_vv_checkpoint

   !> UV file store: `file` or the automatic name, guessing from `guess` or the saved file
   !>
   !> @param[in]  files      UV checkpoint settings
   !> @param[in]  automatic  automatic file name of the model
   !> @param[out] store      file store
   subroutine uv_checkpoint(files, automatic, store)
      !> UV checkpoint settings
      type(rism_checkpoint_files_type), intent(in) :: files
      !> Automatic file name
      character(len=*), intent(in) :: automatic
      !> File store
      type(moz_file_checkpoint_type), intent(out) :: store
      character(len=:), allocatable :: file

      file = automatic
      if (allocated(files%file)) file = trim(files%file)
      if (allocated(files%guess)) then
         call new_file_checkpoint(store, file, [character(len=len_trim(files%guess)) :: trim(files%guess)], &
            & files%interval)
      else
         call new_file_checkpoint(store, file, interval=files%interval)
      end if
   end subroutine uv_checkpoint

   !> Automatic VV checkpoint file, vv_rism_solvent_<id>_<theory>_<closure>.mbin
   !>
   !> - In the directory named by MOIST_RISM_VV_OUTDIR when set, else the working directory
   !>
   !> @param[in] solvent_id  solvent id; 0 for a custom solvent
   !> @param[in] theory      bulk theory name
   !> @param[in] closure     closure name
   function vv_checkpoint_name(solvent_id, theory, closure) result(file)
      !> Solvent id
      integer, intent(in) :: solvent_id
      !> Bulk theory name
      character(len=*), intent(in) :: theory
      !> Closure name
      character(len=*), intent(in) :: closure
      !> File path
      character(len=:), allocatable :: file
      character(len=:), allocatable :: dir
      character(len=16) :: id

      write (id, "(i0)") solvent_id
      file = "vv_rism_solvent_"//trim(id)//"_"//to_lower(trim(adjustl(theory)))//"_"// &
         & to_lower(trim(adjustl(closure)))//".mbin"
      dir = get_env(vv_outdir_env)
      if (len(dir) > 0) file = dir//"/"//file
   end function vv_checkpoint_name

   !> Require the requested LJ/charge suppliers to have actually covered this side
   !>
   !> @param[in] potential     updated, parameter-only potential
   !> @param[out] error        missing supplier or inconsistent solute charge
   !> @param[in] total_charge  optional solute total charge, e
   subroutine check_sites(potential, error, total_charge)
      !> Updated, parameter-only potential
      type(moz_potential_type), intent(in) :: potential
      !> Missing supplier or inconsistent solute charge
      type(error_type), allocatable, intent(out) :: error
      !> Optional solute total charge, e
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
end module moist_app_rism_setup
