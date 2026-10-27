!> RISM input files and the input constructors of the UV models, without solving VV or UV
module test_moz_input
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_context, only: moist_context_type, new_context
   use moist_model_type, only: solvation_model_type
   use moist_model_moz_input, only: model_moz_input_type, moz_checkpoint_files_type, new_vv_checkpoint, &
      & new_uv_checkpoint, vv_checkpoint_name
   use moist_model_moz_1d_input, only: model_moz_1d_input_type, new_moz_1d_model
   use moist_model_moz_3d_input, only: model_moz_3d_input_type, new_moz_3d_model
   use moist_model_moz_1d_type, only: model_moz_1d_type
   use moist_model_moz_3d_type, only: model_moz_3d_type
   use moist_model_moz_potential, only: lj_mixing_geometric, potential_sites_type
   use moist_model_moz_checkpoint_file, only: moz_file_checkpoint_type
   use moist_model_moz_solvent_vv, only: solvent_vv_parameters_type
   use moist_model_moz_solvent_water, only: new_spce_water_geometry
   use moist_math_grid_3d_factory, only: new_spatial_grid
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, &
      & partition_becke, partition_ssf, partition_pvoronoi
   use moist_solvent_system, only: solvent_system_type
   use test_moz_fixtures, only: water_id, new_water_system
   use test_helpers, only: check_moist_error
   implicit none(type, external)
   private
   public :: collect_moz_input

contains

   !> Register input-boundary regression cases
   subroutine collect_moz_input(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("native_schemas", check_native_schemas), &
         & new_unittest("strict_keys", check_strict_keys), &
         & new_unittest("grid_groups", check_grid_groups), &
         & new_unittest("molecular_partitions", check_molecular_partitions), &
         & new_unittest("invalid_grids", check_invalid_grids), &
         & new_unittest("water_presets", check_water_presets), &
         & new_unittest("construction", check_construction), &
         & new_unittest("solvent_settings", check_solvent_settings), &
         & new_unittest("lj_tables", check_lj_tables), &
         & new_unittest("charge_choices", check_charge_choices), &
         & new_unittest("ec_settings", check_ec_settings), &
         & new_unittest("checkpoint_settings", check_checkpoint_settings), &
         & new_unittest("checkpoint_names", check_checkpoint_names)]
   end subroutine collect_moz_input

   !> Read independent TOML, retain separate UV/VV controls and copy optional arrays
   subroutine check_native_schemas(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(model_moz_input_type), allocatable :: input, copy
      integer :: unit, stat, dim
      character(len=*), parameter :: path = "moist-test-moz-input-schemas.toml"

      do dim = 1, 2
         call new_input(input, dim)
         open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write (unit, "(a)") 'closure = "HNC"', "[potential]", "q = [-1.0]", "[convergence]", "tolerance = 2e-8", &
            & "[solvent]", 'closure = "KH"', 'theory = "DRISM"', &
            & "[solvent.drism]", "smear = 0.71", &
            & "[solvent.convergence]", "tolerance = 3e-11"
         close (unit)
         call input%read_file(path, err, strict=.true.)
         call remove_input(path)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, allocated(input%potential%q))
         if (allocated(error)) return
         select type (input)
         type is (model_moz_1d_input_type)
            call check(error, input%uv%closure == "HNC" .and. input%uv%solvent%closure == "KH" &
               & .and. input%uv%tolerance == 2.0e-8_wp .and. input%uv%solvent%tolerance == 3.0e-11_wp &
               & .and. input%uv%solvent%theory == "DRISM" &
               & .and. input%uv%solvent%drism_smear == 0.71_wp &
               & .and. input%grid%kind == "radial")
         type is (model_moz_3d_input_type)
            call check(error, input%uv%closure == "HNC" .and. input%uv%solvent%closure == "KH" &
               & .and. input%uv%tolerance == 2.0e-8_wp .and. input%uv%solvent%tolerance == 3.0e-11_wp &
               & .and. input%uv%solvent%theory == "DRISM" &
               & .and. input%uv%solvent%drism_smear == 0.71_wp &
               & .and. input%grid%kind == "cartesian")
         end select
         if (allocated(error)) return
         allocate (copy, source=input)
         copy%potential%q(1) = 1.0_wp
         call check(error, input%potential%q(1), -1.0_wp)
         if (allocated(error)) return
         deallocate (copy)
         call input%write_file(path, err)
         if (.not. allocated(err)) call input%read_file(path, err, strict=.true.)
         call remove_input(path)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, input%potential%q(1), -1.0_wp)
         if (allocated(error)) return
         call input%init_defaults()
         call check(error,.not. allocated(input%potential%q))
         if (allocated(error)) return
      end do
   end subroutine check_native_schemas

   !> Reject unknown keys, removed keys, empty groups, quoted dotted keys and 3D-only input in 1D
   !>
   !> - Removed: the `[system]` table, the charge `environment` and `lj_fallback`
   !> - EC fit settings only with `charges = "ec"`
   subroutine check_strict_keys(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_1d_input_type) :: input
      integer :: unit, stat, i
      character(len=*), parameter :: path = "moist-test-moz-input-strict.toml"
      character(len=40), parameter :: bad(10) = [character(len=40) :: &
         & 'correction = "PC+"', 'finite_size = "huenenberger"', "unknown = 2", &
         & "[convergence_typo]", '"convergence.tolerance" = 1e-8', "[grid]", "[system]", "[potential]", &
         & "[potential.ec]", "[potential]"]
      character(len=40), parameter :: second(10) = [character(len=40) :: &
         & "", "", "", "", "", "nx = 16", 'solvent = "water"', 'environment = "gas"', "range_scale = 2.0", &
         & 'lj_fallback = "uff"']

      do i = 1, size(bad)
         call input%init_defaults()
         open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write (unit, "(a)") trim(bad(i))
         if (len_trim(second(i)) > 0) write (unit, "(a)") trim(second(i))
         close (unit)
         call input%read_file(path, err, strict=.true.)
         call remove_input(path)
         call check(error, allocated(err), more="accepted: "//trim(bad(i))//" "//trim(second(i)))
         if (allocated(error)) return
         call check(error, index(err%message, "Unknown") > 0 .or. index(err%message, "Literal dots") > 0, &
            & more=err%message)
         if (allocated(error)) return
      end do
   end subroutine check_strict_keys

   !> Keep independent grid settings through TOML input, JSON output and default restoration
   !>
   !> @param[out] error Test failure
   subroutine check_grid_groups(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(model_moz_input_type), allocatable :: input
      integer :: unit, stat, dim
      character(len=9), parameter :: kinds(2) = ["radial   ", "cartesian"]
      character(len=*), parameter :: path = "moist-test-moz-input-grids.toml"
      character(len=*), parameter :: json_path = "moist-test-moz-input-grids.json"

      do dim = 1, 2
         call new_input(input, dim)
         open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write (unit, "(a)", iostat=stat) "[potential]", 'charges = "none"', &
            & "[grid]", 'kind = "'//trim(kinds(dim))//'"', &
            & "[grid.radial]", "npts = 64", "dr = 0.2", &
            & "[grid.cartesian]", "nx = 12", "ny = 14", "nz = 16", "dr = 0.9", "margin = 3.0", &
            & "[grid.molecular]", "dr = 1.1", 'partition = "SSF"', &
            & "[grid.molecular.radial]", "npts = 24", "rmax = 8.0", "handymod_exponent = 0.7", &
            & "[grid.molecular.angular]", "degree = 31"
         call check(error, stat, 0)
         if (allocated(error)) return
         close (unit, iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call input%read_file(path, err, strict=.true.)
         call remove_input(path)
         if (.not. allocated(err)) call input%write_file(json_path, err)
         if (.not. allocated(err)) call input%read_file(json_path, err, strict=.true.)
         call remove_input(json_path)
         call check_moist_error(error, err)
         if (allocated(error)) return
         call check(error, input%grid%kind == trim(kinds(dim)) &
            & .and. input%grid%radial%npts == 64 .and. input%grid%radial%dr == 0.2_wp &
            & .and. input%grid%cartesian%nx == 12 .and. input%grid%cartesian%ny == 14 &
            & .and. input%grid%cartesian%nz == 16 .and. input%grid%cartesian%dr == 0.9_wp &
            & .and. input%grid%cartesian%margin == 3.0_wp &
            & .and. input%grid%molecular%radial%npts == 24 .and. input%grid%molecular%radial%rmax == 8.0_wp &
            & .and. input%grid%molecular%radial%handymod_exponent == 0.7_wp &
            & .and. input%grid%molecular%angular%degree == 31 .and. input%grid%molecular%dr == 1.1_wp &
            & .and. input%grid%molecular%partition == "ssf", "Grid groups lost their independent settings")
         if (allocated(error)) return
         call input%init_defaults()
         call check(error, input%grid%kind == trim(kinds(dim)) &
            & .and. input%grid%radial%dr == 0.015_wp .and. input%grid%cartesian%dr == 0.6_wp &
            & .and. input%grid%molecular%dr == 0.6_wp .and. input%grid%molecular%partition == "becke" &
            & .and. input%grid%molecular%radial%npts == 50 .and. input%grid%molecular%radial%rmax == 10.0_wp &
            & .and. input%grid%molecular%radial%handymod_exponent == 2.0_wp &
            & .and. input%grid%molecular%angular%degree == 17, &
            & "Grid defaults retained earlier input values")
         if (allocated(error)) return
      end do
   end subroutine check_grid_groups

   !> Build all molecular partitions from input and retain their controls through a round trip
   !>
   !> @param[out] error Test failure
   subroutine check_molecular_partitions(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_3d_input_type) :: input
      class(moist_math_grid_3d_type), allocatable :: grid, wider
      type(structure_type) :: mol
      integer :: unit, stat, variant
      integer, parameter :: schemes(3) = [partition_becke, partition_ssf, partition_pvoronoi]
      character(len=8), parameter :: names(3) = ["BECKE   ", "SsF     ", "pvoronoi"]
      character(len=*), parameter :: path = "moist-test-moz-input-partitions.toml"

      call new(mol, [1, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.5_wp, 0.0_wp, 0.0_wp], [3, 2]))
      do variant = 1, size(schemes)
         open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         write (unit, "(a)", iostat=stat) "[potential]", 'charges = "none"', &
            & "[grid]", 'kind = "Molecular"', "[grid.molecular]", &
            & 'partition = "'//trim(names(variant))//'"', "becke_k = 4", "ssf_a = 0.5", "power_width = 0.8", &
            & "weight_pruning_cutoff = 2e-14", "partition_pruning_cutoff = 1e-9", "dr = 1.2", "kbuffer = 1.5", &
            & "nufft_tolerance = 2e-6", "[grid.molecular.radial]", "npts = 8", "rmax = 3.0", &
            & "handymod_exponent = 1.0", "[grid.molecular.angular]", "degree = 17"
         call check(error, stat, 0)
         if (allocated(error)) return
         close (unit, iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call input%read_file(path, err, strict=.true.)
         if (.not. allocated(err)) call input%write_file(path, err)
         if (.not. allocated(err)) call input%read_file(path, err, strict=.true.)
         call remove_input(path)
         if (.not. allocated(err)) call new_spatial_grid(input%grid, grid, err)
         if (.not. allocated(err)) call grid%update(mol, err)
         call check_moist_error(error, err)
         if (allocated(error)) return
         select type (grid)
         type is (moist_math_grid_3d_molecular_type)
            call check(error, grid%get_partition(), schemes(variant))
            if (allocated(error)) return
            select case (schemes(variant))
            case (partition_becke)
               call check(error, grid%get_becke_k(), 4)
            case (partition_ssf)
               call check(error, input%grid%molecular%ssf_a, 0.5_wp)
            case (partition_pvoronoi)
               call check(error, grid%get_power_width(), 0.8_wp)
            case default
               call check(error, .false., "Unknown partition in molecular input test")
            end select
            if (allocated(error)) return
            call check(error, grid%ngrid > 0 .and. all(grid%nrad_per_atom == 8) .and. all(grid%nang_per_atom == 110) &
               & .and. grid%get_dr() == 1.2_wp .and. grid%get_kbuffer() == 1.5_wp &
               & .and. grid%get_nufft_tol() == 2.0e-6_wp .and. grid%get_weight_threshold() == 2.0e-14_wp &
               & .and. grid%get_pruning_threshold() == 1.0e-9_wp, "Molecular grid did not use the input settings")
         class default
            call check(error, .false., "Molecular input constructed another grid kind")
         end select
         if (allocated(error)) return
         if (schemes(variant) == partition_ssf) then
            input%grid%molecular%ssf_a = 1.0_wp
            call new_spatial_grid(input%grid, wider, err)
            if (.not. allocated(err)) call wider%update(mol, err)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, wider%ngrid > grid%ngrid, "Increasing ssf_a did not expand the retained grid")
            if (allocated(error)) return
            deallocate (wider)
         end if
         deallocate (grid)
      end do
   end subroutine check_molecular_partitions

   !> Reject invalid partition controls, grid sizes, thresholds and obsolete flat keys
   !>
   !> @param[out] error Test failure
   subroutine check_invalid_grids(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_3d_input_type) :: input
      integer :: unit, stat, i
      character(len=*), parameter :: path = "moist-test-moz-input-invalid-grids.toml"
      character(len=40), parameter :: bad(25) = [character(len=40) :: &
         & 'partition = "unknown"', "npts = 0", "degree = -1", "rmax = -1.0", &
         & "handymod_exponent = 0.0", "dr = 0.0", "kbuffer = -1.0", "nufft_tolerance = 0.0", &
         & "becke_k = 0", "ssf_a = 0.0", "ssf_a = 1.1", "power_width = 0.0", &
         & "weight_pruning_cutoff = -1.0", "partition_pruning_cutoff = -1.0", "partition_pruning_cutoff = 1.0", &
         & "nx = 1", "npts = 1", "dr = 0.6", "nrad = 50", &
         & "nrad = 50", "angular_degree = 17", "rmax = 10.0", "handymod_exponent = 2.0", &
         & "weight_threshold = 1e-14", "pruning_threshold = 0.0"]

      do i = 1, size(bad)
         open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         select case (i)
         case (2, 4, 5)
            write (unit, "(a)", iostat=stat) "[grid.molecular.radial]", trim(bad(i))
         case (3)
            write (unit, "(a)", iostat=stat) "[grid.molecular.angular]", trim(bad(i))
         case (16)
            write (unit, "(a)", iostat=stat) "[grid.cartesian]", trim(bad(i))
         case (17)
            write (unit, "(a)", iostat=stat) "[grid.radial]", trim(bad(i))
         case (18:19)
            write (unit, "(a)", iostat=stat) "[grid]", trim(bad(i))
         case default
            write (unit, "(a)", iostat=stat) "[grid.molecular]", trim(bad(i))
         end select
         call check(error, stat, 0)
         if (allocated(error)) return
         close (unit, iostat=stat)
         call check(error, stat, 0)
         if (allocated(error)) return
         call input%read_file(path, err, strict=.true.)
         call remove_input(path)
         call check(error, allocated(err), more="accepted invalid grid setting: "//trim(bad(i)))
         if (allocated(error)) return
         call check(error, index(err%message, "grid") > 0, more=err%message)
         if (allocated(error)) return
      end do
   end subroutine check_invalid_grids

   !> Read both water LJ presets with SPC/E charges, preserve them on round trip and build both models
   subroutine check_water_presets(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      class(model_moz_input_type), allocatable :: input
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(structure_type) :: mol
      type(potential_sites_type) :: sites
      class(solvation_model_type), allocatable :: model
      integer :: unit, dim, variant
      character(len=5), parameter :: presets(2) = ["spce ", "cspce"]
      character(len=*), parameter :: path = "moist-test-moz-input-water-presets.toml"

      call new_context(ctx, verbosity=0, nthreads=1)
      call new_water_system(system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new_spce_water_geometry(mol)
      do dim = 1, 2
         do variant = 1, size(presets)
            call new_input(input, dim)
            call check(error, size(input%solvent_potential%lj) == 4 .and. input%solvent_potential%lj(1) == "cspce" &
               & .and. size(input%solvent_potential%charges) == 2 .and. input%solvent_potential%charges(1) == "spce" &
               & .and. input%solvent_potential%charges(2) == "resp-solvent", &
               & "Default solvent must select cSPC/E and SPC/E first")
            if (allocated(error)) return
            open (newunit=unit, file=path, status="replace", action="write")
            write (unit, "(a)") "[potential]", &
               & 'lj = "'//trim(presets(variant))//'"', 'charges = "spce"', &
               & "[solvent.potential]", 'lj = "'//trim(presets(variant))//'"', 'charges = "spce"'
            close (unit)
            call input%read_file(path, err, strict=.true.)
            if (.not. allocated(err)) call input%write_file(path, err)
            if (.not. allocated(err)) call input%read_file(path, err, strict=.true.)
            call remove_input(path)
            call check_moist_error(error, err)
            if (allocated(error)) return
            call check(error, size(input%potential%lj) == 1 .and. input%potential%lj(1) == trim(presets(variant)) &
               & .and. size(input%potential%charges) == 1 .and. input%potential%charges(1) == "spce" &
               & .and. size(input%solvent_potential%lj) == 1 &
               & .and. input%solvent_potential%lj(1) == trim(presets(variant)) &
               & .and. size(input%solvent_potential%charges) == 1 .and. input%solvent_potential%charges(1) == "spce")
            if (allocated(error)) return
            call small_grid(input)
            call input%new_model(ctx, system, model, err)
            if (.not. allocated(err)) call model%update(mol, err)
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
            deallocate (model)
         end do
      end do
   end subroutine check_water_presets

   !> The input overloads build models that wait for update; update checks the solute charge
   !>
   !> - The input is left as given; the constructor validates a copy
   !> - A charge mismatch or a wrong charge count leaves the model not updated
   !> - Inputs without restored defaults are refused
   subroutine check_construction(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_1d_input_type) :: input_1d
      type(model_moz_3d_input_type) :: input_3d, raw
      type(model_moz_1d_type) :: model_1d
      type(model_moz_3d_type) :: model_3d
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(structure_type) :: mol

      call new_context(ctx, verbosity=0, nthreads=1)
      call new_water_system(system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new(mol, [17], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=-1.0_wp)

      call input_1d%init_defaults()
      input_1d%potential%q = [-1.0_wp]
      input_1d%solvent_potential%mixing = "Geometric"
      call small_grid(input_1d)
      call new_moz_1d_model(model_1d, ctx, input_1d, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error,.not. model_1d%is_updated() .and. model_1d%has_solvent() .and. .not. model_1d%solvent_solved() &
         & .and. model_1d%rgrid%npts == 32 .and. model_1d%potential%mixing() == lj_mixing_geometric &
         & .and. input_1d%solvent_potential%mixing == "Geometric")
      if (allocated(error)) return
      call model_1d%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model_1d%is_updated() .and. model_1d%atom_count() == 1)
      if (allocated(error)) return
      mol%charge = 0.0_wp
      call model_1d%update(mol, err)
      call check(error, allocated(err), more="an inconsistent molecular charge was accepted in 1D")
      if (allocated(error)) return
      call check(error, index(err%message, "Charges sum to -1.000000 e") > 0, more=err%message)
      if (allocated(error)) return
      call check(error,.not. model_1d%is_updated() .and. .not. model_1d%potential%is_updated() &
         & .and. model_1d%atom_count() == 0, more="a rejected charge sum left the 1D model updated")
      if (allocated(error)) return
      mol%charge = -1.0_wp

      ! One charge too many for the structure
      input_1d%potential%q = [-1.0_wp, 0.0_wp]
      call new_moz_1d_model(model_1d, ctx, input_1d, system, err)
      if (.not. allocated(err)) call model_1d%update(mol, err)
      call check(error, allocated(err), more="two charges for one atom were accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "2 values") > 0 .and. .not. model_1d%is_updated(), more=err%message)
      if (allocated(error)) return

      call input_3d%init_defaults()
      input_3d%potential%q = [-1.0_wp]
      input_3d%solvent_potential%mixing = "geometric"
      call small_grid(input_3d)
      call new_moz_3d_model(model_3d, ctx, input_3d, system, err)
      if (.not. allocated(err)) call model_3d%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model_3d%is_updated() .and. model_3d%atom_count() == 1 .and. model_3d%has_solvent() &
         & .and. .not. model_3d%solvent_solved() .and. model_3d%grid%ngrid == 4096 &
         & .and. model_3d%potential%mixing() == lj_mixing_geometric)
      if (allocated(error)) return

      ! Charges summing to -1 for a neutral structure
      mol%charge = 0.0_wp
      call model_3d%update(mol, err)
      call check(error, allocated(err), more="an inconsistent molecular charge was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "Charges sum to -1.000000 e") > 0, more=err%message)
      if (allocated(error)) return

      call new_moz_3d_model(model_3d, ctx, raw, system, err)
      call check(error, allocated(err), more="an input without defaults was accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "init_defaults") > 0, more=err%message)
   end subroutine check_construction

   !> `[solvent]` name and state point build the solvent system; a host solvent system takes precedence
   !>
   !> - Without a name the host supplies the solvent
   !> - A failed polymorphic construction leaves the model unallocated
   subroutine check_solvent_settings(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_1d_input_type) :: input_1d
      type(model_moz_3d_input_type) :: input_3d
      type(model_moz_1d_type) :: model_1d
      class(solvation_model_type), allocatable :: model
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system, water
      integer :: unit, stat
      character(len=*), parameter :: path = "moist-test-moz-input-solvent.toml"

      call new_context(ctx, verbosity=0, nthreads=1)
      call new_water_system(water, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call input_1d%init_defaults()
      call input_1d%new_system(system, err)
      call check(error, allocated(err), more="a solvent system without a name was built")
      if (allocated(error)) return
      call check(error, index(err%message, "solvent.name") > 0, more=err%message)
      if (allocated(error)) return

      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)") "[potential]", "q = [-1.0]", "[solvent]", 'name = "methanol"', "temperature = 310.0", &
         & "pressure = 2.0e5", 'closure = "KH"', "[solvent.potential]", 'lj = "uff"', 'charges = "resp-solvent"'
      close (unit)
      call input_1d%read_file(path, err, strict=.true.)
      call remove_input(path)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(input_1d%solvent_state%name) .and. allocated(input_1d%solvent_state%temperature) &
         & .and. allocated(input_1d%solvent_state%pressure), more="solvent name and state point not read")
      if (allocated(error)) return
      call check(error, input_1d%solvent_state%name == "methanol" .and. input_1d%solvent_state%temperature == 310.0_wp &
         & .and. input_1d%solvent_state%pressure == 2.0e5_wp .and. input_1d%uv%solvent%closure == "KH")
      if (allocated(error)) return
      call input_1d%new_system(system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, system%temperature == 310.0_wp .and. system%pressure_si == 2.0e5_wp)
      if (allocated(error)) return

      ! Methanol of the file has six sites; the water of the host has three
      call small_grid(input_1d)
      call new_moz_1d_model(model_1d, ctx, input_1d, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model_1d%solvent_nsite(), 6)
      if (allocated(error)) return
      ! The default solvent lists fall back to GAFF2 and RESP charges for methanol
      call input_1d%solvent_potential%init_defaults()
      call new_moz_1d_model(model_1d, ctx, input_1d, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model_1d%solvent_nsite(), 6)
      if (allocated(error)) return
      call new_moz_1d_model(model_1d, ctx, input_1d, water, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, model_1d%solvent_nsite(), 3)
      if (allocated(error)) return

      ! Invalid state points
      input_1d%solvent_state%temperature = -1.0_wp
      call input_1d%validate(err)
      call check(error, allocated(err), more="a negative temperature was accepted")
      if (allocated(error)) return
      deallocate (input_1d%solvent_state%temperature)
      input_1d%solvent_state%name = " "
      call input_1d%new_system(system, err)
      call check(error, allocated(err), more="a blank solvent name was accepted")
      if (allocated(error)) return

      ! Polymorphic construction from the file's solvent
      call input_3d%init_defaults()
      input_3d%solvent_state%name = "water"
      input_3d%potential%q = [-1.0_wp]
      call small_grid(input_3d)
      call input_3d%new_model(ctx, model, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(model), more="no model was built")
      if (allocated(error)) return
      deallocate (model)
      input_3d%solvent_state%name = "no-such-solvent"
      call input_3d%new_model(ctx, model, err)
      call check(error, allocated(err) .and. .not. allocated(model), more="an unknown solvent built a model")
   end subroutine check_solvent_settings

   !> LJ tables in order of precedence; a bare name is a one-entry list
   !>
   !> - Each atom takes the first table covering it: GAFF2 types methane, UFF fills argon
   !> - Empty lists, unknown or repeated tables and custom LJ beside another table are refused
   subroutine check_lj_tables(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_1d_input_type) :: input
      type(model_moz_1d_type) :: model
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(structure_type) :: mol
      type(potential_sites_type) :: sites
      integer :: unit, stat, i
      character(len=*), parameter :: path = "moist-test-moz-input-lj.toml"
      character(len=16), parameter :: refused(2, 4) = reshape([character(len=16) :: &
         & "uff", "uff", "custom", "uff", "gaff", "", "uff", "lj"], [2, 4])

      call new_context(ctx, verbosity=0, nthreads=1)
      call new_water_system(system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      call input%init_defaults()
      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)") "[potential]", 'lj = ["GAFF2", " uff"]', 'charges = "none"', &
         & "[solvent.potential]", 'lj = "cspce"'
      close (unit)
      call input%read_file(path, err, strict=.true.)
      if (.not. allocated(err)) call input%write_file(path, err)
      if (.not. allocated(err)) call input%read_file(path, err, strict=.true.)
      call remove_input(path)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(input%potential%lj) == 2 .and. size(input%solvent_potential%lj) == 1, &
         & more="LJ lists not kept")
      if (allocated(error)) return
      call check(error, input%potential%lj(1) == "gaff2" .and. input%potential%lj(2) == "uff" &
         & .and. input%solvent_potential%lj(1) == "cspce")
      if (allocated(error)) return

      ! Methane with a distant argon atom, in bohr
      call new(mol, [6, 1, 1, 1, 1, 18], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.19_wp, 1.19_wp, 1.19_wp, &
         & -1.19_wp, -1.19_wp, 1.19_wp, -1.19_wp, 1.19_wp, -1.19_wp, 1.19_wp, -1.19_wp, -1.19_wp, &
         & 12.0_wp, 0.0_wp, 0.0_wp], [3, 6]))
      call small_grid(input)
      call new_moz_1d_model(model, ctx, input, system, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (.not. allocated(err)) call model%potential%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(index(sites%pair%label(:5), "/GAFF") > 0) .and. index(sites%pair%label(6), "/UFF") > 0, &
         & more="GAFF2 then UFF: "//sites%pair%label(1)//" "//sites%pair%label(6))
      if (allocated(error)) return
      input%potential%lj = [character(len=16) :: "uff", "gaff2"]
      call new_moz_1d_model(model, ctx, input, system, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (.not. allocated(err)) call model%potential%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(index(sites%pair%label, "/UFF") > 0), more="UFF first: "//sites%pair%label(1))
      if (allocated(error)) return

      do i = 1, size(refused, 2)
         call input%init_defaults()
         input%potential%charges = [character(len=24) :: "none"]
         input%potential%lj = pack(refused(:, i), refused(:, i) /= "")
         if (refused(2, i) == "lj") input%potential%lj = [character(len=16) ::]
         call input%validate(err)
         call check(error, allocated(err), more="LJ tables accepted: "//refused(1, i)//refused(2, i))
         if (allocated(error)) return
      end do
   end subroutine check_lj_tables

   !> Charge choices: custom, none, spce, ec and the full table names; solute and solvent sides
   subroutine check_charge_choices(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_3d_input_type) :: input
      type(model_moz_1d_input_type) :: input_1d
      integer :: i, unit, stat
      character(len=*), parameter :: path = "moist-test-moz-input-charges.toml"
      character(len=24), parameter :: solute_refused(4) = [character(len=24) :: &
         & "fixed", "cspce", "resp-solvent", "resp"]
      character(len=24), parameter :: solvent_refused(4) = [character(len=24) :: &
         & "ec", "cspce", "hirshfeld", "mbis-vacuum"]
      character(len=24), parameter :: lists_refused(2, 3) = reshape([character(len=24) :: &
         & "spce", "spce", "spce", "ec", "spce", "custom"], [2, 3])

      call input%init_defaults()
      input%potential%charges = [character(len=24) :: "none"]
      input%solvent_potential%charges = [character(len=24) :: "Hirshfeld-Conductor"]
      call input%validate(err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, input%solvent_potential%charges(1) == "hirshfeld-conductor", &
         & more="table charge name: "//trim(input%solvent_potential%charges(1)))
      if (allocated(error)) return

      do i = 1, size(solute_refused)
         call input%init_defaults()
         input%potential%charges = [character(len=24) :: solute_refused(i)]
         call input%validate(err)
         call check(error, allocated(err), more="solute charges accepted: "//trim(solute_refused(i)))
         if (allocated(error)) return
      end do
      do i = 1, size(solvent_refused)
         call input%init_defaults()
         input%potential%charges = [character(len=24) :: "none"]
         input%solvent_potential%charges = [character(len=24) :: solvent_refused(i)]
         call input%validate(err)
         call check(error, allocated(err), more="solvent charges accepted: "//trim(solvent_refused(i)))
         if (allocated(error)) return
      end do

      ! Charge lists keep their order; whole-side choices stand alone
      call input%init_defaults()
      input%potential%charges = [character(len=24) :: "none"]
      input%solvent_potential%charges = [character(len=24) :: " SPCE", "Hirshfeld-Solvent"]
      call input%validate(err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, input%solvent_potential%charges(1) == "spce" &
         & .and. input%solvent_potential%charges(2) == "hirshfeld-solvent")
      if (allocated(error)) return
      do i = 1, size(lists_refused, 2)
         call input%init_defaults()
         input%potential%charges = [character(len=24) :: "none"]
         input%solvent_potential%charges = lists_refused(:, i)
         call input%validate(err)
         call check(error, allocated(err), more="solvent charges accepted: "//lists_refused(1, i)//lists_refused(2, i))
         if (allocated(error)) return
      end do
      call input%init_defaults()
      input%potential%charges = [character(len=24) :: "none", "spce"]
      call input%validate(err)
      call check(error, allocated(err), more="neutral solute charges beside SPC/E were accepted")
      if (allocated(error)) return
      call check(error, index(err%message, "only charge choice") > 0, more=err%message)
      if (allocated(error)) return
      call input%init_defaults()
      input%potential%charges = [character(len=24) ::]
      call input%validate(err)
      call check(error, allocated(err), more="an empty charge list was accepted")
      if (allocated(error)) return

      ! Custom charges need q, and q needs custom charges
      call input%init_defaults()
      call input%validate(err)
      call check(error, allocated(err), more="custom charges without q were accepted")
      if (allocated(error)) return
      input%potential%q = [0.0_wp]
      input%potential%charges = [character(len=24) :: "spce"]
      call input%validate(err)
      call check(error, allocated(err), more="q without custom charges was accepted")
      if (allocated(error)) return

      ! EC fits on a volume grid only
      call input_1d%init_defaults()
      input_1d%potential%charges = [character(len=24) :: "ec"]
      call input_1d%validate(err)
      call check(error, allocated(err), more="EC charges were accepted in 1D")
      if (allocated(error)) return
      call check(error, index(err%message, "volume grid") > 0, more=err%message)
      if (allocated(error)) return

      ! EC fit settings on the solvent side reach validation, which names the solute side
      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)") "[potential]", 'charges = "none"', "[solvent.potential]", 'charges = "ec"', &
         & "[solvent.potential.ec]", "range_scale = 2.0"
      close (unit)
      call input%read_file(path, err, strict=.true.)
      call remove_input(path)
      call check(error, allocated(err), more="EC charges were accepted on the solvent side")
      if (allocated(error)) return
      call check(error, index(err%message, "solute potential") > 0, more=err%message)
   end subroutine check_charge_choices

   !> EC fit settings are read with `charges = "ec"` and reach the EC term of the model
   subroutine check_ec_settings(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_3d_input_type) :: input
      type(model_moz_3d_type) :: model
      type(moist_context_type), target :: ctx
      type(solvent_system_type) :: system
      type(structure_type) :: mol
      integer :: unit, stat
      character(len=*), parameter :: path = "moist-test-moz-input-ec.toml"

      call new_context(ctx, verbosity=0, nthreads=1)
      call new_water_system(system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new(mol, [17], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=-1.0_wp)
      call input%init_defaults()
      call check(error, input%potential%ec%exclusion_radius == 1.0_wp &
         & .and. input%potential%ec%restraint_strength == 0.01_wp, more="EC defaults differ from the EC term")
      if (allocated(error)) return
      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)") "[potential]", 'charges = "ec"', "[potential.ec]", "exclusion_radius = 1.5", &
         & "range_scale = 4.0", "long_range_bias = 0.5", "restraint_strength = -1.0", "charge_penalty = 0.2"
      close (unit)
      call input%read_file(path, err, strict=.true.)
      call remove_input(path)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, input%potential%ec%exclusion_radius == 1.5_wp .and. input%potential%ec%range_scale == 4.0_wp &
         & .and. input%potential%ec%long_range_bias == 0.5_wp .and. input%potential%ec%restraint_strength == -1.0_wp &
         & .and. input%potential%ec%charge_penalty == 0.2_wp, more="EC fit settings not read")
      if (allocated(error)) return

      ! The negative restraint reaches the EC term, which refuses it at update
      call small_grid(input)
      call new_moz_3d_model(model, ctx, input, system, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call model%update(mol, err)
      call check(error, allocated(err), more="a negative EC restraint reached no EC term")
      if (allocated(error)) return
      call check(error, index(err%message, "EC fit requires positive") > 0, more=err%message)
      if (allocated(error)) return

      input%potential%ec%restraint_strength = 0.02_wp
      call new_moz_3d_model(model, ctx, input, system, err)
      if (.not. allocated(err)) call model%update(mol, err)
      call check_moist_error(error, err)
   end subroutine check_ec_settings

   !> Read the checkpoint group and reject invalid names and intervals
   subroutine check_checkpoint_settings(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(model_moz_3d_input_type) :: input
      integer :: unit, stat
      character(len=*), parameter :: path = "moist-test-moz-input-checkpoint.toml"

      call input%init_defaults()
      call check(error, input%checkpoint%enabled .and. input%checkpoint%vv%interval == 0 &
         & .and. input%checkpoint%uv%interval == 0 .and. .not. allocated(input%checkpoint%vv%file) &
         & .and. .not. allocated(input%checkpoint%uv%guess) .and. .not. allocated(input%checkpoint%vv%directory), &
         & more="unexpected checkpoint defaults")
      if (allocated(error)) return
      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      call check(error, stat, 0)
      if (allocated(error)) return
      write (unit, "(a)") "[potential]", "q = [-1.0]", "[checkpoint]", "enabled = false", "[checkpoint.vv]", &
         & 'file = "water.mbin"', "interval = 3", 'directory = "cache"', &
         & "[checkpoint.uv]", 'guess = "start.mbin"', "interval = 2"
      close (unit)
      call input%read_file(path, err, strict=.true.)
      call remove_input(path)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error,.not. input%checkpoint%enabled .and. input%checkpoint%vv%interval == 3 &
         & .and. input%checkpoint%uv%interval == 2, more="checkpoint switch or intervals not read")
      if (allocated(error)) return
      call check(error, allocated(input%checkpoint%vv%file) .and. allocated(input%checkpoint%uv%guess) &
         & .and. allocated(input%checkpoint%vv%directory) .and. .not. allocated(input%checkpoint%uv%file), &
         & more="checkpoint files not read")
      if (allocated(error)) return
      call check(error, input%checkpoint%vv%file == "water.mbin" .and. input%checkpoint%uv%guess == "start.mbin" &
         & .and. input%checkpoint%vv%directory == "cache")
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
      deallocate (input%checkpoint%uv%file)
      input%checkpoint%uv%directory = ""
      call input%validate(err)
      call check(error, allocated(err), more="a blank checkpoint directory was accepted")
   end subroutine check_checkpoint_settings

   !> Automatic names in the checkpoint directory; explicit files and guesses as given
   subroutine check_checkpoint_names(error)
      type(error_type), allocatable, intent(out) :: error
      type(moz_checkpoint_files_type) :: files
      type(solvent_vv_parameters_type) :: params
      type(moz_file_checkpoint_type) :: store
      character(len=*), parameter :: stem = "vv_rism_solvent_175_drism_"

      call check(error, vv_checkpoint_name(water_id, "DRISM", " PSE3") == stem//"pse3.mbin", &
         & more=vv_checkpoint_name(water_id, "DRISM", " PSE3"))
      if (allocated(error)) return

      ! Automatic VV name in a directory, then the other closures' files in the same directory
      call files%init_defaults()
      files%directory = "cache"
      call params%init_defaults()
      params%theory = "DRISM"
      params%closure = "KH"
      call new_vv_checkpoint(files, water_id, params, store)
      call check(error, store%file == "cache/"//stem//"kh.mbin", more=store%file)
      if (allocated(error)) return
      call check(error, size(store%guesses) == 5, more="expected the own file and four other closures")
      if (allocated(error)) return
      call check(error, trim(store%guesses(1)) == store%file &
         & .and. trim(store%guesses(2)) == "cache/"//stem//"hnc.mbin", more=trim(store%guesses(2)))
      if (allocated(error)) return

      ! Explicit file: no directory, no other closures
      files%file = "water.mbin"
      call new_vv_checkpoint(files, water_id, params, store)
      call check(error, store%file == "water.mbin" .and. size(store%guesses) == 1, more=store%file)
      if (allocated(error)) return

      ! Automatic UV name in a directory with a trailing separator
      call files%init_defaults()
      files%directory = "out/"
      call new_uv_checkpoint(files, "uv_rism_3d.mbin", store)
      call check(error, store%file == "out/uv_rism_3d.mbin" .and. size(store%guesses) == 1, more=store%file)
   end subroutine check_checkpoint_names

   !> Allocate and default one input; 1 selects 1D, 2 selects 3D
   !>
   !> @param[out] input  input with restored defaults
   !> @param[in]  dim    1 for 1D, 2 for 3D
   subroutine new_input(input, dim)
      !> Input with restored defaults
      class(model_moz_input_type), allocatable, intent(out) :: input
      !> 1 for 1D, 2 for 3D
      integer, intent(in) :: dim
      if (dim == 1) then
         allocate (model_moz_1d_input_type :: input)
      else
         allocate (model_moz_3d_input_type :: input)
      end if
      call input%init_defaults()
   end subroutine new_input

   !> Small radial and Cartesian grids for construction tests
   !>
   !> @param[inout] input  input settings
   subroutine small_grid(input)
      !> Input settings
      class(model_moz_input_type), intent(inout) :: input
      input%grid%radial%npts = 32
      input%grid%radial%dr = 1.0_wp
      input%grid%cartesian%nx = 16
      input%grid%cartesian%ny = 16
      input%grid%cartesian%nz = 16
      input%grid%cartesian%dr = 1.0_wp
      input%grid%cartesian%margin = 2.0_wp
   end subroutine small_grid

   !> Remove a fixture owned by this test
   !> @param[in] path Fixture path
   subroutine remove_input(path)
      character(len=*), intent(in) :: path
      integer :: unit, stat
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine remove_input
end module test_moz_input
