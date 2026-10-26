!> Standalone RISM input and native construction, without solving VV or UV
module test_rism_input
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_app_rism_input, only: rism_input_type, new_rism_input
   use moist_app_rism_setup, only: new_rism_model, vv_checkpoint_name
   use moist_context, only: moist_context_type, new_context
   use moist_model_type, only: solvation_model_type
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_1d_type, only: model_moz_1d_type
   use moist_model_moz_3d_type, only: model_moz_3d_type
   use moist_model_moz_potential, only: lj_mixing_geometric, potential_sites_type
   use moist_model_moz_solvent_water, only: new_spce_water_geometry
   use test_helpers, only: check_moist_error
   implicit none(type, external)
   private
   public :: collect_rism_input

contains

   !> Register input-boundary regression cases
   subroutine collect_rism_input(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("native_schemas", check_native_schemas), &
         & new_unittest("strict_keys", check_strict_keys), &
         & new_unittest("water_presets", check_water_presets), &
         & new_unittest("construction", check_construction), &
         & new_unittest("checkpoint_settings", check_checkpoint_settings)]
   end subroutine collect_rism_input

   !> Read independent TOML, retain separate UV/VV controls and copy optional arrays
   subroutine check_native_schemas(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rism_input_type) :: input, copy
      integer :: unit, stat, mode
      character(len=6), parameter :: modes(2) = ["rism1d", "rism3d"]
      character(len=*), parameter :: path = "moist-test-rism-input-schemas.toml"

      do mode = 1, 2
         call new_rism_input(input, modes(mode), err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         open(newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write(unit, "(a)") 'closure = "HNC"', "[system]", 'solvent = "water"', "charge = -1.0", &
            & "temperature = 310.0", "[potential]", "q = [-1.0]", "[convergence]", "tolerance = 2e-8", &
            & "[solvent]", 'closure = "KH"', 'theory = "DRISM"', &
            & "[solvent.drism]", "smear = 0.71", &
            & "[solvent.convergence]", "tolerance = 3e-11"
         close(unit)
         call input%read_file(path, err, strict=.true.)
         call remove_input(path)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, allocated(input%system%charge) .and. allocated(input%potential%q))
         if (allocated(error)) return
         call check(error, input%system%charge, -1.0_wp)
         if (allocated(error)) return
         call check(error, input%system%temperature, 310.0_wp)
         if (allocated(error)) return
         select type (params => input%model_parameters)
         type is (model_moz_1d_parameters_type)
            call check(error, params%closure == "HNC" .and. params%solvent%closure == "KH" &
               & .and. params%tolerance == 2.0e-8_wp .and. params%solvent%tolerance == 3.0e-11_wp &
               & .and. params%solvent%theory == "DRISM" &
               & .and. params%solvent%drism_smear == 0.71_wp &
               & .and. input%grid%kind == "radial")
         type is (model_moz_3d_parameters_type)
            call check(error, params%closure == "HNC" .and. params%solvent%closure == "KH" &
               & .and. params%tolerance == 2.0e-8_wp .and. params%solvent%tolerance == 3.0e-11_wp &
               & .and. params%solvent%theory == "DRISM" &
               & .and. params%solvent%drism_smear == 0.71_wp &
               & .and. input%grid%kind == "cartesian")
         end select
         if (allocated(error)) return
         copy = input
         copy%potential%q(1) = 1.0_wp
         call check(error, input%potential%q(1), -1.0_wp)
         if (allocated(error)) return
         call input%write_file(path, err)
         if (.not. allocated(err)) call input%read_file(path, err, strict=.true.)
         call remove_input(path)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, input%potential%q(1), -1.0_wp)
         if (allocated(error)) return
         call check(error, input%system%charge, -1.0_wp)
         if (allocated(error)) return
         call input%init_defaults()
         call check(error, .not. allocated(input%system%charge) .and. .not. allocated(input%potential%q))
         if (allocated(error)) return
      end do
   end subroutine check_native_schemas

   !> Reject unknown keys, empty groups, quoted dotted keys and 3D-only input in 1D
   subroutine check_strict_keys(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rism_input_type) :: input
      integer :: unit, stat, i
      character(len=*), parameter :: path = "moist-test-rism-input-strict.toml"
      character(len=72), parameter :: bad(6) = [character(len=72) :: &
         & 'correction = "PC+"', 'finite_size = "hunenberger"', "unknown = 2", &
         & "[convergence_typo]", '"convergence.tolerance" = 1e-8', "[grid]"]

      do i = 1, size(bad)
         call new_rism_input(input, "rism1d", err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         open(newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write(unit, "(a)") trim(bad(i))
         if (i == size(bad)) write(unit, "(a)") "nx = 16"
         close(unit)
         call input%read_file(path, err, strict=.true.)
         call remove_input(path)
         call check(error, allocated(err))
         if (allocated(error)) return
         call check(error, index(err%message, "Unknown") > 0 .or. index(err%message, "Literal dots") > 0)
         if (allocated(error)) return
      end do
   end subroutine check_strict_keys

   !> Read both named water presets and preserve their choices on round trip
   subroutine check_water_presets(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rism_input_type) :: input
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      type(potential_sites_type) :: sites
      class(solvation_model_type), allocatable :: model
      integer :: unit, mode, variant
      character(len=6), parameter :: modes(2) = ["rism1d", "rism3d"]
      character(len=5), parameter :: presets(2) = ["spce ", "cspce"]
      character(len=*), parameter :: path = "moist-test-rism-water-presets.toml"

      call new_context(ctx, verbosity=0, nthreads=1)
      call new_spce_water_geometry(mol)
      do mode = 1, size(modes)
         do variant = 1, size(presets)
            call new_rism_input(input, modes(mode), err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, input%solvent_potential%lj == "cspce" &
               & .and. input%solvent_potential%charges == "cspce", "Default water must explicitly select cSPC/E")
            if (allocated(error)) return
            open(newunit=unit, file=path, status="replace", action="write")
            write(unit, "(a)") "[system]", 'solvent = "water"', "[potential]", &
               & 'lj = "'//trim(presets(variant))//'"', 'charges = "'//trim(presets(variant))//'"', &
               & "[solvent.potential]", 'lj = "'//trim(presets(variant))//'"', &
               & 'charges = "'//trim(presets(variant))//'"'
            close(unit)
            call input%read_file(path, err, strict=.true.)
            if (.not. allocated(err)) call input%validate(err)
            if (.not. allocated(err)) call input%write_file(path, err)
            if (.not. allocated(err)) call input%read_file(path, err, strict=.true.)
            call remove_input(path)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, input%potential%lj == trim(presets(variant)) &
               & .and. input%potential%charges == trim(presets(variant)) &
               & .and. input%solvent_potential%lj == trim(presets(variant)) &
               & .and. input%solvent_potential%charges == trim(presets(variant)))
            if (allocated(error)) return
            input%grid%npts = 32
            input%grid%nx = 16
            input%grid%ny = 16
            input%grid%nz = 16
            input%grid%dr = 1.0_wp
            input%grid%margin = 2.0_wp
            call new_rism_model(input, ctx, mol, model, err)
            if (.not. allocated(err)) then
               select type (model)
               type is (model_moz_1d_type)
                  call model%potential%sites(sites, err)
               type is (model_moz_3d_type)
                  call model%potential%sites(sites, err)
               end select
            end if
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, all(sites%q == [-0.8476_wp, 0.4238_wp, 0.4238_wp]) &
               & .and. all(sites%pair%covered))
            if (allocated(error)) return
            if (variant == 1) then
               call check(error, all(sites%pair%epsilon(2:3) == 0.0_wp))
            else
               call check(error, all(sites%pair%epsilon(2:3) > 0.0_wp))
            end if
            if (allocated(error)) return
            deallocate(model)
         end do
      end do
   end subroutine check_water_presets

   !> Factories produce updated polymorphic UV models whose bulk remains unsolved
   subroutine check_construction(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rism_input_type) :: input
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(solvation_model_type), allocatable :: model
      integer :: mode
      character(len=6), parameter :: modes(2) = ["rism1d", "rism3d"]

      call new_context(ctx, verbosity=0, nthreads=1)
      call new(mol, [17], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=-1.0_wp)
      do mode = 1, 2
         call new_rism_input(input, modes(mode), err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         input%system%solvent = "water"
         input%potential%q = [-1.0_wp]
         input%solvent_potential%mixing = "geometric"
         input%grid%npts = 32
         input%grid%nx = 16
         input%grid%ny = 16
         input%grid%nz = 16
         input%grid%dr = 1.0_wp
         input%grid%margin = 2.0_wp
         call new_rism_model(input, ctx, mol, model, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, model%is_updated() .and. model%atom_count() == 1)
         if (allocated(error)) return
         select type (model)
         type is (model_moz_1d_type)
            call check(error, model%has_solvent() .and. .not. model%solvent_solved() &
               & .and. model%rgrid%npts == 32 .and. model%potential%mixing() == lj_mixing_geometric)
         type is (model_moz_3d_type)
            call check(error, model%has_solvent() .and. .not. model%solvent_solved() &
               & .and. model%grid%ngrid == 4096 .and. model%potential%mixing() == lj_mixing_geometric)
         end select
         if (allocated(error)) return
         deallocate(model)
         mol%charge = 0.0_wp
         call new_rism_model(input, ctx, mol, model, err)
         call check(error, allocated(err), more="an inconsistent molecular charge was accepted")
         if (allocated(error)) return
         mol%charge = -1.0_wp
      end do
   end subroutine check_construction

   !> Read the checkpoint group, reject invalid files and intervals, name automatic VV files
   subroutine check_checkpoint_settings(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(rism_input_type) :: input
      character(len=:), allocatable :: name
      integer :: unit, stat
      character(len=*), parameter :: path = "moist-test-rism-input-checkpoint.toml"
      character(len=*), parameter :: vv_name = "vv_rism_solvent_175_drism_pse3.mbin"

      call new_rism_input(input, "rism3d", err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, input%checkpoint%enabled .and. input%checkpoint%vv%interval == 5 &
         & .and. input%checkpoint%uv%interval == 5 .and. .not. allocated(input%checkpoint%vv%file) &
         & .and. .not. allocated(input%checkpoint%uv%guess), more="unexpected checkpoint defaults")
      if (allocated(error)) return
      open(newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write(unit, "(a)") "[potential]", "q = [-1.0]", "[checkpoint]", "enabled = false", "[checkpoint.vv]", &
         & 'file = "water.mbin"', "interval = 0", "[checkpoint.uv]", 'guess = "start.mbin"', "interval = 2"
      close(unit)
      call input%read_file(path, err, strict=.true.)
      call remove_input(path)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. input%checkpoint%enabled .and. input%checkpoint%vv%interval == 0 &
         & .and. input%checkpoint%uv%interval == 2, more="checkpoint switch or intervals not read")
      if (allocated(error)) return
      call check(error, allocated(input%checkpoint%vv%file) .and. allocated(input%checkpoint%uv%guess) &
         & .and. .not. allocated(input%checkpoint%uv%file), more="checkpoint files not read")
      if (allocated(error)) return
      call check(error, input%checkpoint%vv%file == "water.mbin" .and. input%checkpoint%uv%guess == "start.mbin")
      if (allocated(error)) return
      call input%validate(err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      input%checkpoint%uv%interval = -1
      call input%validate(err)
      call check(error, allocated(err), more="a negative checkpoint interval was accepted")
      if (allocated(error)) return
      input%checkpoint%uv%interval = 1
      input%checkpoint%uv%file = " "
      call input%validate(err)
      call check(error, allocated(err), more="a blank checkpoint file was accepted")
      if (allocated(error)) return

      ! MOIST_RISM_VV_OUTDIR may prefix a directory; the base name is fixed
      name = vv_checkpoint_name(175, "DRISM", " PSE3")
      call check(error, len(name) >= len(vv_name), more="automatic VV name too short: "//name)
      if (allocated(error)) return
      call check(error, name(len(name) - len(vv_name) + 1:) == vv_name, more="unexpected automatic VV name: "//name)
   end subroutine check_checkpoint_settings

   !> Remove a fixture owned by this test
   !> @param[in] path Fixture path
   subroutine remove_input(path)
      character(len=*), intent(in) :: path
      integer :: unit, stat
      open(newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close(unit, status="delete")
   end subroutine remove_input
end module test_rism_input
