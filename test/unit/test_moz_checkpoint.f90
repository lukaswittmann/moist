!> Tests of the MOZ checkpoints: states, the file and memory stores and warm starts: moz_checkpoint
!>
!> - File store round trip, candidate order, refusals and failed saves on small synthetic
!>   states; a file of another format version is skipped
!> - Saves that a load would reject, and malformed stored layouts, each refused in a test of
!>   its own; a refused save keeps the previous state
!> - The long-range tail of a header is known, zero or unknown, and matters only for
!>   interpolated points
!> - A candidate that covers no point is skipped for the next one; stored layouts that need
!>   no refusal; a point grid without atom-centered data restarts on itself
!> - Memory store: copies share the slot of the caller; the VV solvent reloads its
!>   solution from the slot, and the 1D and 3D UV models keep theirs across a geometry update
!> - Point positions of atom-centered grids: compatibility, guesses for the points
!>   retained at both structures, file round trip, and the 3D UV header of a molecular grid
!> - cSPC/E water VV and water-in-water 1D and 3D UV solves with a memory store;
!>   a solve that iterates saves its final field, so a reload without a new save
!>   was accepted without iterations
!> - Warm starts on other grids: VV and 1D UV states on finer or longer radial grids; 3D UV
!>   states on a larger box, a finer lattice and around a rotated solute, and from a
!>   molecular grid onto Cartesian and molecular grids
!> - VV and 1D UV solves that stop on a closure switch succeed only on the requested closure
!> - `moz_blowup`: the live-grid and layout calls agree; a Cartesian grid never updated is refused
!> - Each test runs as a selected test
module test_moz_checkpoint
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, moist_error => error_type, fatal_error
   use mctc_io, only: structure_type, new
   use mctc_io_codata2018, only: Boltzmann_constant, atomic_unit_of_energy
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_context, only: moist_context_type, new_context
   use moist_model_type, only: solvation_model_type
   use moist_channels_coupling, only: coupling_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, new_molecular_point_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_generator_lebedev_type, &
      & moist_math_grid_angular_request_type, moist_math_grid_angular_type
   use moist_model_moz_checkpoint, only: moz_checkpoint_type, moz_checkpoint_state_type, &
      & checkpoint_kind_vv, checkpoint_kind_uv1d, checkpoint_kind_uv3d, tail_known, tail_zero, tail_unknown, &
      & moz_checkpoint_guess
   use moist_model_moz_potential_kernel_ng, only: ng_real_tail
   use moist_model_moz_blowup, only: moz_blowup
   use moist_math_grid_layout, only: moist_math_grid_layout_type
   use moist_math_grid_transfer, only: grid_transfer_options_type, grid_transfer_report_type
   use moist_utils_binio, only: binary_file
   use moist_model_moz_checkpoint_file, only: moz_file_checkpoint_type, new_file_checkpoint
   use moist_model_moz_checkpoint_memory, only: moz_memory_checkpoint_type, new_memory_checkpoint
   use moist_utils_file, only: create_new_file, delete_file
   use moist_model_moz_solvent_vv, only: solvent_vv_type, solvent_vv_parameters_type, solvent_vv_thermodynamics_type
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_potential, only: lj_cspce, coulomb_spce
   use test_helpers, only: check_moist_error, expect_moist_error, get_qc_handymod_recipe, get_test_water
   use test_moz_fixtures, only: new_cspce_water_solvent
   implicit none(type, external)
   private
   public :: collect_moz_checkpoint

   !> Built-in Lebedev points in reverse order: the same counts, other point positions
   type, extends(moist_math_grid_angular_generator_lebedev_type) :: permuted_angular_type
   contains
      procedure :: select => permuted_angular_select
   end type permuted_angular_type

   !> Saved states of a memory store, shared by its copies
   type :: store_log_type
      !> States in order of saving
      type(moz_checkpoint_state_type), allocatable :: saved(:)
   end type store_log_type

   !> Store keeping states in memory; load offers `preset`
   type, extends(moz_checkpoint_type) :: memory_store_type
      !> Shared log of the saves
      type(store_log_type), pointer :: log => null()
      !> State offered by load; none without fields
      type(moz_checkpoint_state_type) :: preset
   contains
      procedure :: save => memory_save
      procedure :: load => memory_load
      procedure :: describe => memory_describe
   end type memory_store_type

   !> Files of the file-store tests
   character(len=*), parameter :: file_a = "moist-test-checkpoint-a.mbin", file_b = "moist-test-checkpoint-b.mbin", &
      & file_c = "moist-test-checkpoint-c.mbin", file_d = "moist-test-checkpoint-d.mbin", &
      & file_trip = "moist-test-checkpoint-trip.mbin", file_save = "moist-test-checkpoint-save.mbin", &
      & file_points = "moist-test-checkpoint-points.mbin", file_version = "moist-test-checkpoint-version.mbin"
   !> Files of the candidate test of moz_checkpoint_guess, apart from file_candidates' files
   character(len=*), parameter :: file_blowup_a = "moist-test-checkpoint-blowup-a.mbin", &
      & file_blowup_b = "moist-test-checkpoint-blowup-b.mbin", file_blowup_c = "moist-test-checkpoint-blowup-c.mbin"
   !> Files of the refused-save tests
   character(len=*), parameter :: file_no_origin = "moist-test-checkpoint-no-origin.mbin", &
      & file_field_shape = "moist-test-checkpoint-field-shape.mbin", file_keep = "moist-test-checkpoint-keep.mbin"
   !> Point positions of the synthetic molecular state, among 24 unpruned points
   integer, parameter :: stored_points(7) = [1, 2, 5, 7, 8, 13, 20]

contains

   !> Register tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_moz_checkpoint(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("file_round_trip", check_file_round_trip), &
         & new_unittest("file_candidates", check_file_candidates), &
         & new_unittest("file_failed_save", check_file_failed_save), &
         & new_unittest("file_save_no_origin", check_file_save_no_origin, should_fail=.true.), &
         & new_unittest("file_save_field_shape", check_file_save_field_shape, should_fail=.true.), &
         & new_unittest("memory_save_no_nodes", check_memory_save_no_nodes, should_fail=.true.), &
         & new_unittest("memory_save_rows", check_memory_save_rows, should_fail=.true.), &
         & new_unittest("refused_save_keeps", check_refused_save_keeps), &
         & new_unittest("file_version", check_file_version), &
         & new_unittest("tail_status", check_tail_status), &
         & new_unittest("blowup_candidates", check_blowup_candidates), &
         & new_unittest("blowup_tail", check_blowup_tail), &
         & new_unittest("tail_values", check_tail_values), &
         & new_unittest("blowup_layout", check_blowup_layout), &
         & new_unittest("layout_weights_size", check_layout_weights_size, should_fail=.true.), &
         & new_unittest("layout_kmodes_axes", check_layout_kmodes_axes, should_fail=.true.), &
         & new_unittest("layout_box_dims", check_layout_box_dims, should_fail=.true.), &
         & new_unittest("layout_box_corner", check_layout_box_corner, should_fail=.true.), &
         & new_unittest("layout_box_overflow", check_layout_box_overflow, should_fail=.true.), &
         & new_unittest("layout_corner_nan", check_layout_corner_nan, should_fail=.true.), &
         & new_unittest("layout_radial_nodes", check_layout_radial_nodes, should_fail=.true.), &
         & new_unittest("layout_radial_order", check_layout_radial_order, should_fail=.true.), &
         & new_unittest("layout_points_nan", check_layout_points_nan, should_fail=.true.), &
         & new_unittest("layout_no_fields", check_layout_no_fields, should_fail=.true.), &
         & new_unittest("points_restart", check_points_restart), &
         & new_unittest("blowup_overloads", check_blowup_overloads), &
         & new_unittest("blowup_unbuilt_box", check_blowup_unbuilt_box, should_fail=.true.), &
         & new_unittest("ramp_switch_vv", check_ramp_switch_vv), &
         & new_unittest("ramp_switch_uv1d", check_ramp_switch_uv1d), &
         & new_unittest("vv_reload", check_vv_reload), &
         & new_unittest("vv_intermediate", check_vv_intermediate), &
         & new_unittest("vv_ramp", check_vv_ramp), &
         & new_unittest("vv_unfitting", check_vv_unfitting), &
         & new_unittest("uv1d_reload", check_uv1d_reload), &
         & new_unittest("uv3d_reload", check_uv3d_reload), &
         & new_unittest("points_compatible", check_points_compatible), &
         & new_unittest("points_angular", check_points_angular), &
         & new_unittest("points_guess", check_points_guess), &
         & new_unittest("file_points", check_file_points), &
         & new_unittest("uv3d_molecular", check_uv3d_molecular), &
         & new_unittest("memory_store", check_memory_store), &
         & new_unittest("memory_update", check_memory_update), &
         & new_unittest("memory_vv", check_memory_vv), &
         & new_unittest("memory_uv1d", check_memory_uv1d), &
         & new_unittest("transfer_vv", check_transfer_vv), &
         & new_unittest("transfer_uv1d", check_transfer_uv1d), &
         & new_unittest("transfer_uv3d_box", check_transfer_uv3d_box), &
         & new_unittest("transfer_uv3d_finer", check_transfer_uv3d_finer), &
         & new_unittest("transfer_uv3d_rotated", check_transfer_uv3d_rotated), &
         & new_unittest("transfer_molecular", check_transfer_molecular), &
         & new_unittest("transfer_cartesian_molecular", check_transfer_cartesian_molecular), &
         & new_unittest("transfer_uv3d_points", check_transfer_uv3d_points)]
   end subroutine collect_moz_checkpoint

   !> Save a state to a file and read it back unchanged
   !>
   !> - Every dataset survives, including the optional provenance
   !> - An empty file name saves nothing; describe names the saved file
   subroutine check_file_round_trip(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err, warning
      type(moz_file_checkpoint_type) :: store, quiet
      type(moz_checkpoint_state_type) :: state, loaded
      character(len=:), allocatable :: source
      logical :: found

      call synthetic_state(state)
      call new_file_checkpoint(store, file_trip, interval=3)
      call check(error, store%interval == 3 .and. store%describe() == file_trip, more="store settings lost")
      if (allocated(error)) return
      call store%save(state, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      store%header = header_of(state)
      call store%load(loaded, found, source, warning)
      call delete_file(file_trip, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, found .and. .not. allocated(warning) .and. source == file_trip, more="saved state not found")
      if (allocated(error)) return
      call check(error, loaded%kind == state%kind .and. all(loaded%grid%shape == state%grid%shape) &
         & .and. loaded%grid%spacing == state%grid%spacing .and. all(loaded%grid%origin == state%grid%origin) &
         & .and. all(loaded%solvent_sites == state%solvent_sites) .and. all(loaded%solute_sites == state%solute_sites) &
         & .and. all(loaded%solute_xyz == state%solute_xyz) .and. loaded%grid%kind == state%grid%kind, &
         & more="discretization or sites changed")
      if (allocated(error)) return
      call check(error, loaded%tail == state%tail .and. loaded%alpha == state%alpha &
         & .and. all(loaded%solute_charges == state%solute_charges) &
         & .and. all(loaded%solvent_charges == state%solvent_charges), more="long-range tail changed")
      if (allocated(error)) return
      call check(error, loaded%closure == state%closure .and. loaded%theory == state%theory &
         & .and. loaded%temperature == state%temperature .and. loaded%residual == state%residual &
         & .and. loaded%lambda == state%lambda .and. loaded%iterations == state%iterations &
         & .and. (loaded%converged .eqv. state%converged), more="provenance changed")
      if (allocated(error)) return
      call check(error, all(loaded%field_shape == [8, 2]) .and. all(loaded%tr == state%tr), more="t(r) changed")
      if (allocated(error)) return

      call new_file_checkpoint(quiet, "")
      call quiet%save(state, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, len(quiet%describe()) == 0, more="an empty file name describes a location")
   end subroutine check_file_round_trip

   !> Load from candidate files in order
   !>
   !> - Missing files are skipped silently; no candidate leaves nothing found
   !> - Unreadable files, states of another solute, and fields that do not fit their grid
   !>   are skipped with a warning
   !> - The first fitting file is loaded, also as an initial guess
   subroutine check_file_candidates(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err, warning
      type(moist_context_type), target :: ctx
      type(moz_file_checkpoint_type) :: store, writer
      type(moz_checkpoint_state_type) :: state, other, malformed, loaded
      real(wp), allocatable :: guess(:, :)
      character(len=:), allocatable :: source
      character(len=len(file_a)) :: files(5)
      logical :: found
      integer :: unit

      call new_context(ctx, nthreads=0, verbosity=0)
      call synthetic_state(state)
      other = state
      other%solute_sites = [9]
      malformed = state
      malformed%tr = malformed%tr(:7, :)
      call new_file_checkpoint(writer, file_b)
      call writer%save(other, err)
      ! Stores refuse such a state, so it is written as a bare container
      if (.not. allocated(err)) call write_bare_state(file_d, malformed, err)
      if (.not. allocated(err)) then
         call new_file_checkpoint(writer, file_c)
         call writer%save(state, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      open (newunit=unit, file=file_a, status="replace", action="write")
      write (unit, "(a)") "not a container"
      close (unit)

      files = [character(len=len(file_a)) :: "moist-test-checkpoint-0.mbin", file_a, file_b, file_d, file_c]
      call new_file_checkpoint(store, file_c, files)
      store%header = header_of(state)
      store%header%field_shape = [8, 2]
      call store%load(loaded, found, source, warning)
      call check(error, found .and. source == file_c, more="the fitting last candidate was not loaded")
      if (allocated(error)) return
      call check(error, allocated(warning), more="unreadable and unfitting candidates were skipped silently")
      if (allocated(error)) return
      call check(error, index(warning%message, file_a) > 0 .and. index(warning%message, file_b) > 0 .and. &
         & index(warning%message, file_d) > 0 .and. index(warning%message, file_c) == 0, &
         & more="candidate warnings do not identify exactly the rejected files: "//warning%message)
      if (allocated(error)) return
      call check(error, all(loaded%tr == state%tr), more="loaded the wrong fields")
      if (allocated(error)) return

      store%header = header_of(state)
      call store%load_guess(ctx, 1.0_wp, spread([0.0_wp, 0.0_wp], 1, 8), 4, guess)
      call delete_file(file_a, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call delete_file(file_b, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call delete_file(file_d, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, allocated(guess), more="the malformed candidate hid the fitting one from the guess")
      if (allocated(error)) return
      call check(error, all(guess == state%tr), more="the guess holds the wrong fields")
      if (allocated(error)) return

      store%header%solvent_sites = [8, 1, 1]
      call store%load(loaded, found, source, warning)
      call check(error, .not. found .and. allocated(warning), more="a state of other sites was loaded")
      if (allocated(error)) return

      call delete_file(file_c, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      store%header = header_of(state)
      call store%load(loaded, found, source, warning)
      call check(error, .not. found .and. .not. allocated(warning), more="missing files were not skipped silently")
   end subroutine check_file_candidates

   !> Saves replace the file whole and keep it on failure
   !>
   !> - A temporary file held by another save is left untouched
   !> - A later save replaces the existing file
   !> - An incomplete state is refused and the previous file stays loadable
   subroutine check_file_failed_save(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err, warning
      type(moz_file_checkpoint_type) :: store
      type(moz_checkpoint_state_type) :: state, incomplete, loaded
      character(len=:), allocatable :: source, other
      logical :: found, held
      integer :: nbytes

      call synthetic_state(state)
      call new_file_checkpoint(store, file_save)
      store%header = header_of(state)
      call create_new_file(file_save, other, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call store%save(state, err)
      if (.not. allocated(err)) then
         state%residual = 2.5e-10_wp
         call store%save(state, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) then
         call delete_file(file_save, cleanup_error)
         if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
         if (allocated(error)) return
         if (allocated(other)) call delete_file(other, cleanup_error)
         if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
         if (allocated(error)) return
         return
      end if
      inquire (file=other, exist=held, size=nbytes)
      call delete_file(other, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return

      incomplete = state
      deallocate (incomplete%tr)
      incomplete%residual = 1.0_wp
      call store%save(incomplete, err)
      call store%load(loaded, found, source, warning)
      call delete_file(file_save, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, held .and. nbytes == 0, more="the temporary file of another save was changed")
      if (allocated(error)) return
      call check(error, allocated(err), more="an incomplete state was saved")
      if (allocated(error)) return
      call check(error, found .and. .not. allocated(warning), more="the failed save lost the previous file")
      if (allocated(error)) return
      call check(error, loaded%residual == state%residual, more="the later save did not replace the file")
   end subroutine check_file_failed_save

   !> A file store refuses a Cartesian state without its corner
   !>
   !> @param[out] error  the refusal of the save
   subroutine check_file_save_no_origin(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err
      type(moz_file_checkpoint_type) :: store
      type(moz_checkpoint_state_type) :: state

      call synthetic_state(state)
      deallocate (state%grid%origin)
      call new_file_checkpoint(store, file_no_origin)
      call store%save(state, err)
      call delete_file(file_no_origin, cleanup_error)
      if (allocated(cleanup_error)) return
      call expect_moist_error(error, err, "box is not described")
   end subroutine check_file_save_no_origin

   !> A file store refuses a state whose stored field shape differs from its t(r)
   !>
   !> @param[out] error  the refusal of the save
   subroutine check_file_save_field_shape(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err
      type(moz_file_checkpoint_type) :: store
      type(moz_checkpoint_state_type) :: state

      call synthetic_state(state)
      state%field_shape = [8, 3]
      call new_file_checkpoint(store, file_field_shape)
      call store%save(state, err)
      call delete_file(file_field_shape, cleanup_error)
      if (allocated(cleanup_error)) return
      call expect_moist_error(error, err, "field shape differs")
   end subroutine check_file_save_field_shape

   !> A memory store refuses a radial state without its nodes
   !>
   !> @param[out] error  the refusal of the save
   subroutine check_memory_save_no_nodes(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moz_checkpoint_state_type), target :: slot
      type(moz_memory_checkpoint_type) :: store
      type(moz_checkpoint_state_type) :: state

      call synthetic_state(state)
      state%kind = checkpoint_kind_uv1d
      state%grid%kind = "radial"
      state%grid%shape = [8]
      deallocate (state%grid%origin)
      call new_memory_checkpoint(store, slot)
      call store%save(state, err)
      call expect_moist_error(error, err, "radial nodes are not described")
   end subroutine check_memory_save_no_nodes

   !> A memory store refuses a t(r) with fewer rows than the points of its layout
   !>
   !> @param[out] error  the refusal of the save
   subroutine check_memory_save_rows(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moz_checkpoint_state_type), target :: slot
      type(moz_memory_checkpoint_type) :: store
      type(moz_checkpoint_state_type) :: state

      call synthetic_state(state)
      state%tr = state%tr(:7, :)
      call new_memory_checkpoint(store, slot)
      call store%save(state, err)
      call expect_moist_error(error, err, "fields do not match the grid")
   end subroutine check_memory_save_rows

   !> A refused save keeps the previous state of either store loadable
   !>
   !> @param[out] error  test failure
   subroutine check_refused_save_keeps(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err, warning
      type(moz_checkpoint_state_type), target :: slot
      type(moz_memory_checkpoint_type) :: memory
      type(moz_file_checkpoint_type) :: file
      type(moz_checkpoint_state_type) :: state, broken, loaded
      character(len=:), allocatable :: source
      logical :: found

      call synthetic_state(state)
      broken = state
      deallocate (broken%grid%origin)
      broken%residual = 1.0_wp

      call new_memory_checkpoint(memory, slot)
      memory%header = header_of(state)
      call memory%save(state, err)
      call check_moist_error(error, err, "memory save")
      if (allocated(error)) return
      call memory%save(broken, err)
      call check(error, allocated(err), more="a memory store saved a box without its corner")
      if (allocated(error)) return
      call memory%load(loaded, found, source, warning)
      call check(error, found .and. .not. allocated(warning), more="the refused save lost the memory state")
      if (allocated(error)) return
      call check(error, loaded%residual == state%residual, more="the memory slot changed")
      if (allocated(error)) return

      call new_file_checkpoint(file, file_keep)
      file%header = header_of(state)
      call file%save(state, err)
      if (allocated(err)) then
         call delete_file(file_keep, cleanup_error)
         if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
         if (allocated(error)) return
         call check_moist_error(error, err, "file save")
         return
      end if
      call file%save(broken, err)
      call file%load(loaded, found, source, warning)
      call delete_file(file_keep, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, allocated(err), more="a file store saved a box without its corner")
      if (allocated(error)) return
      call check(error, found .and. .not. allocated(warning), more="the refused save lost the file")
      if (allocated(error)) return
      call check(error, loaded%residual == state%residual, more="the file changed")
   end subroutine check_refused_save_keeps

   !> A file of another format version is skipped with a warning
   !>
   !> @param[out] error  test failure
   subroutine check_file_version(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err, warning
      type(moz_file_checkpoint_type) :: store
      type(moz_checkpoint_state_type) :: state, loaded
      type(binary_file) :: binf
      character(len=:), allocatable :: source
      logical :: found

      call binf%open(file_version, status="new", action="w", error=err)
      if (.not. allocated(err)) call binf%add("/checkpoint/version", 1, error=err)
      if (.not. allocated(err)) call binf%close(error=err)
      call check_moist_error(error, err)
      if (allocated(error)) then
         call delete_file(file_version, cleanup_error)
         if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
         if (allocated(error)) return
         return
      end if
      call synthetic_state(state)
      call new_file_checkpoint(store, file_version)
      store%header = header_of(state)
      call store%load(loaded, found, source, warning)
      call delete_file(file_version, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, .not. found .and. allocated(warning), more="a file of format version 1 was loaded")
      if (allocated(error)) return
      call check(error, index(warning%message, "format version") > 0, more="unexpected warning "//warning%message)
   end subroutine check_file_version

   !> The long-range tail of a header is known, zero or unknown
   !>
   !> - The Ng tail rebuilt on the Cartesian points reproduces the potential: known
   !> - A perturbed potential, or a tail without charges: unknown
   !> - No tail and a vanishing potential: zero; a header without points: unknown
   !>
   !> @param[out] error  test failure
   subroutine check_tail_status(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(memory_store_type) :: store
      type(store_log_type), target :: log
      type(moz_checkpoint_state_type) :: state
      real(wp), allocatable :: xyz(:, :), ur_lr(:, :), f(:)
      real(wp), parameter :: alpha = 0.8_wp, q_u(1) = [-1.0_wp], q_v(2) = [-0.8476_wp, 0.4238_wp]
      integer :: v

      call synthetic_state(state)
      call new_memory_store(store, log)
      store%header = header_of(state)
      store%header%grid = state%grid
      store%header%solute_xyz = state%solute_xyz
      call store%header%grid%points(xyz)
      call check(error, size(xyz, 2) == 8, more="Cartesian layout without its points")
      if (allocated(error)) return
      allocate (f(size(xyz, 2)), ur_lr(size(xyz, 2), 2))
      call ng_real_tail(norm2(xyz - spread(state%solute_xyz(:, 1), 2, size(xyz, 2)), dim=1), alpha, f)
      do v = 1, 2
         ur_lr(:, v) = q_v(v)*q_u(1)*f
      end do
      call store%set_tail(ur_lr, alpha, q_u, q_v)
      call check(error, store%header%tail == tail_known .and. store%header%alpha == alpha &
         & .and. all(store%header%solute_charges == q_u) .and. all(store%header%solvent_charges == q_v), &
         & more="reproduced tail not known")
      if (allocated(error)) return
      ur_lr(3, 2) = ur_lr(3, 2) + 1.0e-6_wp
      call store%set_tail(ur_lr, alpha, q_u, q_v)
      call check(error, store%header%tail == tail_unknown .and. .not. allocated(store%header%solute_charges), &
         & more="perturbed tail not unknown")
      if (allocated(error)) return
      call store%set_tail(ur_lr, 0.0_wp)
      call check(error, store%header%tail == tail_unknown, more="tail without charges not unknown")
      if (allocated(error)) return
      ur_lr = 0.0_wp
      call store%set_tail(ur_lr, 0.0_wp)
      call check(error, store%header%tail == tail_zero, more="vanishing tail not zero")
      if (allocated(error)) return
      ur_lr(:, :) = 1.0_wp
      store%header%grid%kind = "radial"
      call store%set_tail(ur_lr, alpha, q_u, q_v)
      call check(error, store%header%tail == tail_unknown, more="tail on a header without points not unknown")
   end subroutine check_tail_status

   !> A candidate that cannot seed the grid is skipped for the next one
   !>
   !> - A state on a point grid far from the current points, and one whose points lie next
   !>   to the current points but in one plane, pass the site checks but cover no point of
   !>   the Cartesian grid: the first is beyond reach, the second has no stencil; the file
   !>   store moves on to a Cartesian state that does
   !>
   !> @param[out] error  test failure
   subroutine check_blowup_candidates(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moz_file_checkpoint_type) :: store, writer
      type(moz_checkpoint_state_type) :: state, points, flat
      real(wp), allocatable :: guess(:, :), xyz(:, :)
      logical, allocatable :: matched(:)
      character(len=len(file_blowup_a)) :: files(3)
      integer :: i, j

      call new_context(ctx, nthreads=0, verbosity=0)
      call synthetic_state(state)
      points = state
      call state%grid%points(xyz)
      deallocate (points%grid%origin)
      points%grid%kind = "points"
      points%grid%shape = [8]
      points%grid%xyz = xyz + 100.0_wp
      allocate (points%grid%w(8), source=0.125_wp)
      allocate (points%grid%owner(8), source=1)
      ! 36 points in the plane between the two layers of the box, 0.29 bohr from its points
      flat = points
      flat%grid%shape = [36]
      deallocate (flat%grid%xyz, flat%grid%w, flat%grid%owner, flat%tr)
      allocate (flat%grid%xyz(3, 36), flat%tr(36, 2))
      do j = 0, 5
         do i = 0, 5
            flat%grid%xyz(:, 1 + i + 6*j) = [-1.85_wp + 0.5_wp*real(i, wp), -2.85_wp + 0.5_wp*real(j, wp), -2.5_wp]
         end do
      end do
      allocate (flat%grid%w(36), source=0.125_wp)
      allocate (flat%grid%owner(36), source=1)
      flat%tr = 1.0_wp
      call new_file_checkpoint(writer, file_blowup_a)
      call writer%save(points, err)
      if (.not. allocated(err)) then
         call new_file_checkpoint(writer, file_blowup_c)
         call writer%save(flat, err)
      end if
      if (.not. allocated(err)) then
         call new_file_checkpoint(writer, file_blowup_b)
         call writer%save(state, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) return
      files = [character(len=len(file_blowup_a)) :: file_blowup_a, file_blowup_c, file_blowup_b]
      call new_file_checkpoint(store, file_blowup_b, files)
      store%header = header_of(state)
      call store%load_guess(ctx, 1.0_wp, spread([0.0_wp, 0.0_wp], 1, 8), 4, guess, matched)
      call delete_file(file_blowup_a, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call delete_file(file_blowup_b, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call delete_file(file_blowup_c, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, allocated(guess), more="the usable third candidate was not tried")
      if (allocated(error)) return
      call check(error, all(matched) .and. all(guess == state%tr), more="the third candidate was not copied")
      call ctx%delete()
   end subroutine check_blowup_candidates

   !> The long-range tail matters only when the stored field is interpolated
   !>
   !> - Unknown tail, coincident lattice: copied
   !> - Unknown tail or incomplete known data, solute moved against the box: refused
   !> - Complete known tail, solute moved against the box: interpolated
   !>
   !> @param[out] error  test failure
   subroutine check_blowup_tail(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(moz_checkpoint_state_type) :: state, header, shifted
      real(wp), allocatable :: guess(:, :)
      logical, allocatable :: matched(:)
      character(len=:), allocatable :: method, reason
      real(wp) :: ur_lr(8, 2)

      call new_context(ctx, nthreads=0, verbosity=0)
      ur_lr = 0.0_wp
      call synthetic_state(state)
      header = header_of(state)
      ! The solute moved against the box: the lattice no longer coincides
      shifted = header
      shifted%solute_xyz(1, 1) = shifted%solute_xyz(1, 1) + 0.1_wp
      state%tail = tail_unknown
      call moz_checkpoint_guess(state, header, ctx, 1.0_wp, ur_lr, guess, matched, method, reason)
      call check(error, method == "copy" .and. allocated(guess), more="unknown tail blocked the copy: "//reason)
      if (allocated(error)) return
      call moz_checkpoint_guess(state, shifted, ctx, 1.0_wp, ur_lr, guess, matched, method, reason)
      call check(error, .not. allocated(guess) .and. index(reason, "unknown") > 0, &
         & more="an unknown tail was interpolated")
      if (allocated(error)) return
      state%tail = tail_known
      deallocate (state%solute_charges)
      call moz_checkpoint_guess(state, shifted, ctx, 1.0_wp, ur_lr, guess, matched, method, reason)
      call check(error, .not. allocated(guess) .and. index(reason, "incomplete") > 0, &
         & more="a known tail without charges was interpolated")
      if (allocated(error)) return
      call synthetic_state(state)
      call moz_checkpoint_guess(state, shifted, ctx, 1.0_wp, ur_lr, guess, matched, method, reason)
      call check(error, method == "fourier" .and. allocated(guess), more="a complete tail was not interpolated: "//reason)
      call ctx%delete()
   end subroutine check_blowup_tail

   !> Interpolation removes the stored Ng tail and restores it at moved points by value
   !>
   !> @param[out] error  test failure
   subroutine check_tail_values(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(moz_checkpoint_state_type) :: state, header
      real(wp), allocatable :: x(:, :), y(:, :), f(:), current_lr(:, :), guess(:, :), expected(:, :)
      logical, allocatable :: matched(:)
      character(len=:), allocatable :: method, reason
      real(wp) :: beta, shift(3), k
      integer :: v

      call new_context(ctx, nthreads=0, verbosity=0)
      call synthetic_state(state)
      state%grid%shape = [8, 8, 8]
      state%grid%origin = -2.0_wp*[1.0_wp, 1.0_wp, 1.0_wp]
      state%solute_xyz = 0.0_wp*state%solute_xyz
      call state%grid%points(x)
      beta = 1.0_wp/(Boltzmann_constant/atomic_unit_of_energy*state%temperature)
      k = 2.0_wp*atan(1.0_wp)
      allocate (f(size(x, 2)))
      call ng_real_tail(norm2(x, dim=1), state%alpha, f)
      deallocate (state%tr)
      allocate (state%tr(size(x, 2), 2))
      do v = 1, 2
         state%tr(:, v) = 0.3_wp*v + 0.2_wp*v*cos(k*x(1, :)) &
            & + beta*state%solute_charges(1)*state%solvent_charges(v)*f
      end do
      header = header_of(state)
      shift = [0.12_wp, -0.08_wp, 0.05_wp]
      header%solute_xyz(:, 1) = shift
      y = x - spread(shift, 2, size(x, 2))
      call ng_real_tail(norm2(y, dim=1), state%alpha, f)
      allocate (current_lr(size(x, 2), 2), expected(size(x, 2), 2))
      do v = 1, 2
         current_lr(:, v) = state%solute_charges(1)*state%solvent_charges(v)*f
         expected(:, v) = 0.3_wp*v + 0.2_wp*v*cos(k*y(1, :))
      end do
      call moz_checkpoint_guess(state, header, ctx, beta, current_lr, guess, matched, method, reason)
      call check(error, allocated(guess) .and. method == "fourier", more="tail transfer refused: "//reason)
      if (allocated(error)) return
      call check(error, all(matched), more="the moved points should all lie inside the source box")
      if (allocated(error)) return
      call check(error, maxval(abs(guess - expected)) < 1.0e-7_wp, &
         & more="the interpolated short-range field contains a missing or interpolated Ng tail")
      call ctx%delete()
   end subroutine check_tail_values

   !> A point grid without atom-centered data restarts on itself
   !>
   !> - The stored points, in another order, are copied onto the current points, also
   !>   with an unknown tail
   !> - Points moved by more than the copy tolerance take no value
   !>
   !> @param[out] error  test failure
   subroutine check_points_restart(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(moz_checkpoint_state_type) :: state, points, header
      real(wp), allocatable :: guess(:, :), xyz(:, :)
      logical, allocatable :: matched(:)
      character(len=:), allocatable :: method, reason
      integer, parameter :: order(8) = [5, 2, 8, 1, 7, 3, 6, 4]
      real(wp) :: ur_lr(8, 2)

      call new_context(ctx, nthreads=0, verbosity=0)
      ur_lr = 0.0_wp
      call synthetic_state(state)
      call state%grid%points(xyz)
      points = state
      deallocate (points%grid%origin)
      points%grid%kind = "points"
      points%grid%shape = [8]
      points%grid%xyz = xyz(:, order)
      allocate (points%grid%w(8), source=0.125_wp)
      allocate (points%grid%owner(8), source=1)
      points%tr = state%tr(order, :)
      header = header_of(points)
      header%grid%xyz = xyz
      call moz_checkpoint_guess(points, header, ctx, 1.0_wp, ur_lr, guess, matched, method, reason)
      call check(error, method == "copy" .and. allocated(guess), more="the shuffled points were not copied: "//reason)
      if (allocated(error)) return
      call check(error, all(matched) .and. all(guess == state%tr), more="copied values in the wrong order")
      if (allocated(error)) return
      points%tail = tail_unknown
      call moz_checkpoint_guess(points, header, ctx, 1.0_wp, ur_lr, guess, matched, method, reason)
      call check(error, method == "copy" .and. allocated(guess), more="an unknown tail blocked the copy: "//reason)
      if (allocated(error)) return
      call check(error, all(guess == state%tr), more="copied values with an unknown tail")
      if (allocated(error)) return
      header%grid%xyz(1, 3) = header%grid%xyz(1, 3) + 1.0e-6_wp
      call moz_checkpoint_guess(points, header, ctx, 1.0_wp, ur_lr, guess, matched, method, reason)
      call check(error, method == "copy" .and. count(matched) == 7 .and. .not. matched(3), &
         & more="a moved point was copied")
      call ctx%delete()
   end subroutine check_points_restart

   !> The live-grid and layout calls of `moz_blowup` agree
   !>
   !> - Cartesian onto a larger, finer Cartesian box and onto a molecular grid; radial onto
   !>   finer and longer nodes
   !> - Points beyond the old grid keep the caller's values
   !>
   !> @param[out] error  test failure
   subroutine check_blowup_overloads(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: small, large
      type(moist_math_grid_3d_molecular_type) :: molecular
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_radial_type) :: coarse, fine, kspace
      type(moist_math_grid_layout_type) :: layout_old, layout_new
      type(grid_transfer_options_type) :: options
      type(grid_transfer_report_type) :: report
      real(wp), allocatable :: f(:, :), g_grid(:, :), g_layout(:, :), exact(:)
      logical, allocatable :: covered_grid(:), covered_layout(:)

      call get_test_water(mol)
      call new_cartesian_point_grid(small, err, 12, 12, 12, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call small%update(mol, err)
      if (.not. allocated(err)) call new_cartesian_point_grid(large, err, 20, 20, 20, 0.4_wp, margin=0.0_wp)
      if (.not. allocated(err)) call large%update(mol, err)
      if (.not. allocated(err)) call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(molecular, err, recipe=recipe, pruning_threshold=0.1_wp, &
         & dr=1.1_wp)
      if (.not. allocated(err)) call molecular%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      allocate (f(small%ngrid, 2))
      call small%get_layout(layout_old)
      f(:, 1) = periodic_field(small%xyz)
      f(:, 2) = -2.0_wp*f(:, 1)
      options%nthreads = 1
      ! Cartesian onto a larger and finer box
      allocate (g_grid(large%ngrid, 2), g_layout(large%ngrid, 2), source=-3.0_wp)
      call moz_blowup(small, f, large, g_grid, err, covered=covered_grid, options=options, report=report)
      call small%get_layout(layout_old)
      call large%get_layout(layout_new)
      if (.not. allocated(err)) call moz_blowup(layout_old, f, layout_new, g_layout, err, covered=covered_layout, options=options)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(g_grid == g_layout) .and. all(covered_grid .eqv. covered_layout), &
         & more="grid and layout calls differ, Cartesian")
      if (allocated(error)) return
      call check(error, report%ninterpolated > 0 .and. any(.not. covered_grid) .and. &
         & all(pack(g_grid, spread(.not. covered_grid, 2, 2)) == -3.0_wp), more="points beyond the old box")
      if (allocated(error)) return
      exact = periodic_field(large%xyz)
      call check(error, maxval(abs(pack(g_grid(:, 1) - exact, covered_grid))), 0.0_wp, thr=1.0e-8_wp, &
         & more="blowup differs from the analytic periodic field on covered large points")
      if (allocated(error)) return
      call check(error, maxval(abs(pack(g_grid(:, 2) + 2.0_wp*g_grid(:, 1), covered_grid))) < 1.0e-12_wp, &
         & more="blowup mixed the two field columns")
      if (allocated(error)) return
      ! Cartesian onto a molecular grid
      deallocate (g_grid, g_layout)
      allocate (g_grid(molecular%ngrid, 2), g_layout(molecular%ngrid, 2), source=0.0_wp)
      call moz_blowup(small, f, molecular, g_grid, err, covered=covered_grid, options=options)
      call molecular%get_layout(layout_new)
      if (.not. allocated(err)) call moz_blowup(layout_old, f, layout_new, g_layout, err, covered=covered_layout, options=options)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(g_grid == g_layout) .and. all(covered_grid .eqv. covered_layout) .and. any(covered_grid), &
         & more="grid and layout calls differ, molecular")
      if (allocated(error)) return
      exact = periodic_field(molecular%xyz)
      call check(error, maxval(abs(pack(g_grid(:, 1) - exact, covered_grid))) < 0.03_wp, &
         & more="blowup differs from the analytic periodic field on covered molecular points")
      if (allocated(error)) return
      call check(error, maxval(abs(pack(g_grid(:, 2) + 2.0_wp*g_grid(:, 1), covered_grid))) < 1.0e-12_wp, &
         & more="blowup mixed the two field columns")
      if (allocated(error)) return
      ! Radial onto finer and longer nodes
      call new_uniform_radial_pair(coarse, kspace, 64, 0.1_wp, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(fine, kspace, 256, 0.03_wp, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      deallocate (f, g_grid, g_layout)
      allocate (f(64, 1), g_grid(256, 1), g_layout(256, 1), source=0.0_wp)
      f(:, 1) = exp(-coarse%r**2)*cos(coarse%r)
      call moz_blowup(coarse, f, fine, g_grid, err, covered=covered_grid, report=report)
      call coarse%get_layout(layout_old)
      call fine%get_layout(layout_new)
      if (.not. allocated(err)) call moz_blowup(layout_old, f, layout_new, g_layout, err, covered=covered_layout)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(g_grid == g_layout) .and. all(covered_grid .eqv. covered_layout) .and. &
         & any(.not. covered_grid) .and. index(report%method, "spline") > 0, more="grid and layout calls differ, radial")
      if (allocated(error)) return
      exact = exp(-fine%r**2)*cos(fine%r)
      call check(error, maxval(abs(pack(g_grid(:, 1) - exact, covered_grid))) < 1.0e-4_wp, &
         & more="radial blowup differs from the analytic field")
   contains
      !> A smooth field periodic on the source box
      !>
      !> @param[in] xyz  points (3, npoint), bohr
      pure function periodic_field(xyz) result(values)
         !> Evaluation points
         real(wp), intent(in) :: xyz(:, :)
         !> Analytic field at the points
         real(wp) :: values(size(xyz, 2))
         real(wp) :: k

         k = 8.0_wp*atan(1.0_wp)/(real(small%nx, wp)*small%dr)
         values = 1.0_wp + 0.2_wp*cos(k*(xyz(1, :) - layout_old%origin(1))) &
            & + 0.1_wp*sin(k*(xyz(2, :) - layout_old%origin(2)))
      end function periodic_field
   end subroutine check_blowup_overloads

   !> Refused: `moz_blowup` from a Cartesian grid that was never updated
   !>
   !> @param[out] error  the refusal
   subroutine check_blowup_unbuilt_box(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moist_math_grid_3d_cartesian_type) :: unbuilt, box
      real(wp), allocatable :: f(:, :), g(:, :)

      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new_cartesian_point_grid(unbuilt, err, 8, 8, 8, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_cartesian_point_grid(box, err, 8, 8, 8, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call box%update(mol, err)
      if (allocated(err)) return
      allocate (f(8**3, 1), source=1.0_wp)
      allocate (g(box%ngrid, 1), source=0.0_wp)
      call moz_blowup(unbuilt, f, box, g, err)
      if (.not. allocated(err)) return
      if (index(err%message, "box is not described") == 0) return
      call check_moist_error(error, err)
   end subroutine check_blowup_unbuilt_box

   !> Stored layouts that need no refusal
   !>
   !> - Other radial nodes of the same count are splined, not copied
   !> - A state with neither a field shape nor fields seeds no guess
   !> - A known tail whose screening parameter is not finite is incomplete
   !>
   !> @param[out] error  test failure
   subroutine check_blowup_layout(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error

      call radial_case()
      if (allocated(error)) return
      call empty_case()
      if (allocated(error)) return
      call tail_cases()
   contains
      !> A state without field shape and fields seeds no guess
      subroutine empty_case()
         type(moist_context_type), target :: ctx
         type(moz_checkpoint_state_type) :: state, header
         real(wp), allocatable :: guess(:, :), ur_lr(:, :)
         logical, allocatable :: matched(:)
         character(len=:), allocatable :: method, reason

         call molecular_state(state, stored_points)
         header = header_of(state)
         allocate (ur_lr(size(state%tr, 1), size(state%tr, 2)), source=0.0_wp)
         deallocate (state%tr)
         if (allocated(state%field_shape)) deallocate (state%field_shape)
         call new_context(ctx, nthreads=0, verbosity=0)
         call moz_checkpoint_guess(state, header, ctx, 1.0_wp, ur_lr, guess, matched, method, reason)
         call ctx%delete()
         call check(error, .not. allocated(guess) .and. len(method) == 0 .and. len(reason) > 0, &
            & more="a state without fields seeded a guess")
      end subroutine empty_case

      !> Other radial nodes of the same count are no copy
      subroutine radial_case()
         type(moist_context_type), target :: ctx
         type(moz_checkpoint_state_type) :: radial, target
         real(wp), allocatable :: guess(:, :)
         logical, allocatable :: matched(:)
         character(len=:), allocatable :: method, reason
         integer :: i

         call radial_state(radial, target)
         target%grid%r = [(0.25_wp*real(i, wp)**1.5_wp, i=1, 6)]
         call new_context(ctx, nthreads=0, verbosity=0)
         call moz_checkpoint_guess(radial, target, ctx, 1.0_wp, reshape([(0.0_wp, i=1, 6)], [6, 1]), &
            & guess, matched, method, reason)
         call ctx%delete()
         call check(error, method == "spline", more="other radial nodes of the same count copied: "//method)
      end subroutine radial_case

      !> A known tail with a screening parameter that is not finite is incomplete
      subroutine tail_cases()
         type(moist_context_type), target :: ctx
         type(moz_checkpoint_state_type) :: cart, moved
         real(wp), allocatable :: guess(:, :)
         logical, allocatable :: matched(:)
         character(len=:), allocatable :: method, reason
         real(wp) :: ur_lr(8, 2)

         call new_context(ctx, nthreads=0, verbosity=0)
         ur_lr = 0.0_wp
         call synthetic_state(cart)
         moved = header_of(cart)
         moved%solute_xyz(1, 1) = moved%solute_xyz(1, 1) + 0.1_wp
         cart%alpha = ieee_value(1.0_wp, ieee_quiet_nan)
         call moz_checkpoint_guess(cart, moved, ctx, 1.0_wp, ur_lr, guess, matched, method, reason)
         call check(error, .not. allocated(guess) .and. index(reason, "incomplete") > 0, &
            & more="a tail with alpha not finite was interpolated: "//reason)
         call ctx%delete()
      end subroutine tail_cases
   end subroutine check_blowup_layout

   !> VV state on six radial nodes and its header, with the field shape
   !>
   !> @param[out] state   state with fields (6, 1)
   !> @param[out] header  the same nodes without fields
   subroutine radial_state(state, header)
      !> State
      type(moz_checkpoint_state_type), intent(out) :: state
      !> Header
      type(moz_checkpoint_state_type), intent(out) :: header
      integer :: i

      state%kind = checkpoint_kind_vv
      state%solvent_sites = [8, 1, 1]
      state%grid%kind = "radial"
      state%grid%shape = [6]
      state%grid%spacing = 0.0_wp
      state%grid%r = [(0.3_wp*real(i, wp)**1.5_wp, i=1, 6)]
      allocate (state%tr(6, 1), source=1.0_wp)
      state%field_shape = [6, 1]
      header = state
      deallocate (header%tr)
   end subroutine radial_state

   !> A malformed stored layout is refused before its fields are read
   !>
   !> @param[out] error    the refusal
   !> @param[in]  variant  defect: 1 weights, 2 Fourier basis, 3 box dimensions, 4 corner,
   !>                      5 overflowing dimensions, 6 corner not finite, 7 radial nodes,
   !>                      8 node order, 9 coordinates not finite, 10 no fields
   subroutine stored_layout_refused(error, variant)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Defect
      integer, intent(in) :: variant
      type(moist_error), allocatable :: err
      type(moz_checkpoint_state_type) :: state, header
      !> Part of the expected refusal
      character(len=:), allocatable :: expected

      select case (variant)
      case (1, 2, 9, 10)
         call molecular_state(state, stored_points)
         header = header_of(state)
      case (7, 8)
         call radial_state(state, header)
      case default
         call synthetic_state(state)
         header = header_of(state)
      end select
      select case (variant)
      case (1)
         state%grid%w = state%grid%w(2:)
         expected = "weights do not match"
      case (2)
         state%grid%kmodes = [4, 4]
         expected = "Fourier basis is not described"
      case (3)
         state%grid%shape = [2, 2, 3]
         expected = "fields do not match the grid"
      case (4)
         deallocate (state%grid%origin)
         expected = "box is not described"
      case (5)
         ! The product wraps to the eight field rows without the overflow guard.
         state%grid%shape = [2, 2, 1073741826]
         expected = "box is not described"
      case (6)
         state%grid%origin(2) = ieee_value(1.0_wp, ieee_quiet_nan)
         expected = "box is not described"
      case (7)
         deallocate (state%grid%r)
         expected = "radial nodes are not described"
      case (8)
         state%grid%r = [1.0_wp, 0.5_wp, 1.5_wp, 2.0_wp, 2.5_wp, 3.0_wp]
         expected = "radial nodes are not described"
      case (9)
         state%grid%xyz(2, 3) = ieee_value(1.0_wp, ieee_quiet_nan)
         expected = "points are not described"
      case default
         deallocate (state%tr)
         if (allocated(state%field_shape)) deallocate (state%field_shape)
         expected = "no indirect correlation"
      end select
      call state%check_compatible(header, err)
      call expect_moist_error(error, err, expected)
   end subroutine stored_layout_refused

   !> Refused: weights of another size than the points
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_weights_size(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 1)
   end subroutine check_layout_weights_size

   !> Refused: a Fourier basis of two axes
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_kmodes_axes(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 2)
   end subroutine check_layout_kmodes_axes

   !> Refused: a box of other dimensions than its fields
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_box_dims(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 3)
   end subroutine check_layout_box_dims

   !> Refused: a box without its corner
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_box_corner(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 4)
   end subroutine check_layout_box_corner

   !> Refused: box dimensions whose product overflows
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_box_overflow(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 5)
   end subroutine check_layout_box_overflow

   !> Refused: a box corner that is not finite
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_corner_nan(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 6)
   end subroutine check_layout_corner_nan

   !> Refused: a radial state without nodes
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_radial_nodes(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 7)
   end subroutine check_layout_radial_nodes

   !> Refused: radial nodes that do not increase
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_radial_order(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 8)
   end subroutine check_layout_radial_order

   !> Refused: point coordinates that are not finite
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_points_nan(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 9)
   end subroutine check_layout_points_nan

   !> Refused: a state with neither a field shape nor fields
   !>
   !> @param[out] error  the refusal
   subroutine check_layout_no_fields(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      call stored_layout_refused(error, 10)
   end subroutine check_layout_no_fields

   !> A VV solve that stops on a closure switch succeeds only on the requested closure
   !>
   !> - A tolerance far above every residual makes the HNC ramp switch after each step
   !> - One step ends on KH -> PSE2: the solve fails as unconverged
   !> - Three steps end on PSE3 -> HNC: the solve succeeds and saves the HNC residual
   !>
   !> @param[out] error  test failure
   subroutine check_ramp_switch_vv(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(solvent_vv_parameters_type) :: params
      type(solvent_vv_type) :: solvent
      real(wp), parameter :: loose = 1.0e10_wp
      real(wp) :: saved_residual

      call new_context(ctx, nthreads=0, verbosity=0)
      call vv_settings(params)
      params%closure = "HNC"
      params%closure_ramp = .true.
      params%tolerance = loose
      params%maxiter = 1
      call new_memory_store(store, log)
      call solve_vv(ctx, params, store, solvent, err)
      call check(error, allocated(err), more="a VV solve stopped on an intermediate closure succeeded")
      if (allocated(error)) return
      call check(error, index(err%message, "1D VV solve did not converge") > 0, more=err%message)
      if (allocated(error)) return
      call check(error, size(log%saved) == 0, more="a field of an intermediate closure was saved")
      if (allocated(error)) return
      params%maxiter = 3
      call solve_vv(ctx, params, store, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the field after the switch to HNC was not saved")
      if (allocated(error)) return
      call check(error, log%saved(1)%residual > 1.0e-6_wp .and. log%saved(1)%residual < loose &
         & .and. log%saved(1)%closure == "HNC", &
         & more="the saved residual is not that of HNC")
      if (allocated(error)) return
      saved_residual = log%saved(1)%residual
      store%preset = log%saved(1)
      call solve_vv(ctx, params, store, solvent, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the accepted ramp checkpoint iterated")
      if (allocated(error)) return
      call check(error, solvent%residual(), saved_residual, rel=.true., thr=1.0e-8_wp, &
         & more="the saved ramp residual differs from the HNC residual of the reloaded field")
      call ctx%delete()
   end subroutine check_ramp_switch_vv

   !> A 1D UV solve that stops on a closure switch succeeds only on the requested closure
   !>
   !> @param[out] error  test failure
   subroutine check_ramp_switch_uv1d(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(model_moz_1d_parameters_type) :: params
      type(model_moz_1d_type) :: model
      real(wp) :: energy
      real(wp), parameter :: loose = 1.0e10_wp
      real(wp) :: saved_residual

      call new_context(ctx, nthreads=0, verbosity=0)
      call vv_settings(params%solvent)
      params%closure = "HNC"
      params%closure_ramp = .true.
      params%tolerance = loose
      params%maxiter = 1
      call new_memory_store(store, log)
      energy = 0.0_wp
      call new_water_1d(ctx, params, store, model, err)
      if (.not. allocated(err)) call model_energy(model, energy, err)
      call check(error, allocated(err), more="a 1D UV solve stopped on an intermediate closure succeeded")
      if (allocated(error)) return
      call check(error, index(err%message, "1D UV solve did not converge") > 0, more=err%message)
      if (allocated(error)) return
      params%maxiter = 3
      energy = 0.0_wp
      call new_water_1d(ctx, params, store, model, err)
      if (.not. allocated(err)) call model_energy(model, energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the field after the switch to HNC was not saved")
      if (allocated(error)) return
      call check(error, log%saved(1)%residual > 1.0e-6_wp .and. log%saved(1)%residual < loose &
         & .and. log%saved(1)%closure == "HNC", &
         & more="the saved residual is not that of HNC")
      if (allocated(error)) return
      saved_residual = log%saved(1)%residual
      store%preset = log%saved(1)
      params%tolerance = 1.001_wp*saved_residual
      call new_water_1d(ctx, params, store, model, err)
      energy = 0.0_wp
      if (.not. allocated(err)) call model_energy(model, energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the ramp checkpoint failed its looser HNC tolerance")
      if (allocated(error)) return
      params%tolerance = 0.999_wp*saved_residual
      call new_water_1d(ctx, params, store, model, err)
      energy = 0.0_wp
      if (.not. allocated(err)) call model_energy(model, energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 2, more="the ramp checkpoint met a tighter HNC tolerance without iterating")
      if (allocated(error)) return
      call check(error, log%saved(2)%iterations > 0, more="the tighter HNC tolerance took no steps")
      call ctx%delete()
   end subroutine check_ramp_switch_uv1d

   !> Reload a converged VV solution without iterations
   !>
   !> - The cold solve saves one converged state at full coupling
   !> - A new solvent loading it saves nothing and has the same h(r) and mu
   subroutine check_vv_reload(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(solvent_vv_parameters_type) :: params
      type(solvent_vv_type) :: cold, warm
      real(wp) :: mu_cold, mu_warm, dh

      call new_context(ctx, nthreads=0, verbosity=0)
      call vv_settings(params)
      call new_memory_store(store, log)
      call solve_vv(ctx, params, store, cold, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the cold solve did not save exactly its final field")
      if (allocated(error)) return
      associate (saved => log%saved(1))
         call check(error, saved%converged .and. saved%iterations > 0 .and. saved%lambda == 1.0_wp &
            & .and. saved%residual < params%tolerance .and. saved%kind == checkpoint_kind_vv &
            & .and. all(saved%grid%shape == [params%npts]) .and. all(saved%solvent_sites == [8, 1, 1]) &
            & .and. saved%closure == "KH" .and. saved%grid%kind == "radial" .and. size(saved%grid%r) == params%npts, &
            & more="unexpected saved VV state")
      end associate
      if (allocated(error)) return

      store%preset = log%saved(1)
      call solve_vv(ctx, params, store, warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1 .and. warm%is_solved(params), &
         & more="the converged reload iterated or left the solvent unsolved")
      if (allocated(error)) return
      call vv_compare(cold, warm, mu_cold, mu_warm, dh, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, dh < 1.0e-12_wp, more="reloaded h(r) differs")
      if (allocated(error)) return
      call check(error, mu_warm, mu_cold, thr=1.0e-13_wp, more="reloaded chemical potential differs")
   end subroutine check_vv_reload

   !> Save every iteration and restart from an unconverged state
   !>
   !> - Intermediate states are unconverged and numbered by iteration; the last is converged
   !> - Restarting from a late unconverged state converges in fewer iterations to the same mu
   subroutine check_vv_intermediate(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(solvent_vv_parameters_type) :: params
      type(solvent_vv_type) :: cold, warm
      type(moz_checkpoint_state_type) :: restart
      real(wp) :: mu_cold, mu_warm, dh
      integer :: n, i, cold_iterations

      call new_context(ctx, nthreads=0, verbosity=0)
      call vv_settings(params)
      call new_memory_store(store, log)
      store%interval = 1
      call solve_vv(ctx, params, store, cold, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      n = size(log%saved)
      call check(error, n >= 3, more="too few intermediate states were saved")
      if (allocated(error)) return
      do i = 1, n - 1
         call check(error, .not. log%saved(i)%converged .and. log%saved(i)%iterations == i, &
            & more="intermediate state marked converged or misnumbered")
         if (allocated(error)) return
      end do
      call check(error, log%saved(n)%converged .and. log%saved(n)%iterations >= n - 1, &
         & more="final state not converged")
      if (allocated(error)) return
      cold_iterations = log%saved(n)%iterations

      ! A new store and empty log, as for the cold solve: ifort 2021 at -O2 otherwise keeps the
      ! size of a log emptied here across the solve that saves into it
      restart = log%saved(n - 2)
      call new_memory_store(store, log)
      store%interval = 1
      store%preset = restart
      call solve_vv(ctx, params, store, warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      n = size(log%saved)
      call check(error, n >= 1, more="the restart saved nothing")
      if (allocated(error)) return
      call check(error, log%saved(n)%converged .and. log%saved(n)%iterations < cold_iterations, &
         & more="the restart did not converge faster than the cold solve")
      if (allocated(error)) return
      call vv_compare(cold, warm, mu_cold, mu_warm, dh, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, mu_warm, mu_cold, rel=.true., thr=1.0e-6_wp, more="restart reached another solution")
   end subroutine check_vv_intermediate

   !> Closure ramp: saves of the requested closure only, and a guess skips the ramp
   !>
   !> - HNC ramped from KH saves only HNC fields
   !> - A converged guess is accepted without iterations despite the ramp
   subroutine check_vv_ramp(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(solvent_vv_parameters_type) :: params
      type(solvent_vv_type) :: cold, warm
      integer :: n

      call new_context(ctx, nthreads=0, verbosity=0)
      call vv_settings(params)
      params%closure = "HNC"
      call new_memory_store(store, log)
      store%interval = 1
      call solve_vv(ctx, params, store, cold, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      n = size(log%saved)
      call check(error, n >= 1, more="the ramped solve saved nothing")
      if (allocated(error)) return
      ! Each of KH, PSE2 and PSE3 ends with at least one iteration
      call check(error, log%saved(1)%iterations > 3 .and. log%saved(n)%converged, &
         & more="a state below the requested closure was saved")
      if (allocated(error)) return

      store%preset = log%saved(n)
      call solve_vv(ctx, params, store, warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == n .and. warm%is_solved(params), &
         & more="a converged guess ran the closure ramp")
   end subroutine check_vv_ramp

   !> States of other sites or another kind are ignored; the solve starts cold
   subroutine check_vv_unfitting(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(solvent_vv_parameters_type) :: params
      type(solvent_vv_type) :: cold, warm

      call new_context(ctx, nthreads=0, verbosity=0)
      call vv_settings(params)
      call new_memory_store(store, log)
      call solve_vv(ctx, params, store, cold, err)
      call check_moist_error(error, err)
      if (allocated(error)) return

      store%preset = log%saved(1)
      store%preset%solvent_sites = [8, 1]
      call solve_vv(ctx, params, store, warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 2 .and. log%saved(2)%iterations > 0 .and. warm%is_solved(params), &
         & more="a state of other sites was used")
      if (allocated(error)) return

      store%preset = log%saved(1)
      store%preset%kind = checkpoint_kind_uv1d
      call solve_vv(ctx, params, store, warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 3 .and. log%saved(3)%iterations > 0, more="a UV state was used for VV")
   end subroutine check_vv_unfitting

   !> Reload a converged 1D UV solution without iterations
   !>
   !> - The energy saves one converged 1D UV state with the solute sites
   !> - A new model loading it saves nothing and has the same energy
   subroutine check_uv1d_reload(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(model_moz_1d_parameters_type) :: params
      type(model_moz_1d_type), target :: cold, warm
      real(wp) :: e_cold, e_warm

      call new_context(ctx, nthreads=0, verbosity=0)
      params%tolerance = 1.0e-8_wp
      call vv_settings(params%solvent)
      call new_memory_store(store, log)
      call new_water_1d(ctx, params, store, cold, err)
      e_cold = 0.0_wp
      if (.not. allocated(err)) call model_energy(cold, e_cold, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the cold solve did not save exactly its final field")
      if (allocated(error)) return
      associate (saved => log%saved(1))
         call check(error, saved%converged .and. saved%iterations > 0 .and. saved%kind == checkpoint_kind_uv1d &
            & .and. all(saved%solute_sites == [8, 1, 1]) .and. all(saved%grid%shape == [params%solvent%npts]) &
            & .and. saved%grid%kind == "radial" .and. .not. allocated(saved%tail), more="unexpected saved 1D UV state")
      end associate
      if (allocated(error)) return

      store%preset = log%saved(1)
      call new_water_1d(ctx, params, store, warm, err)
      e_warm = 0.0_wp
      if (.not. allocated(err)) call model_energy(warm, e_warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the converged reload iterated")
      if (allocated(error)) return
      call check(error, e_warm, e_cold, thr=1.0e-12_wp, more="reloaded energy differs")
   end subroutine check_uv1d_reload

   !> Reload a converged 3D UV solution, and warm-start a displaced solute
   !>
   !> - The energy saves one converged 3D UV state on the Cartesian grid
   !> - A new model loading it saves nothing and has the same energy
   !> - The same state starts a slightly displaced solute, which converges in
   !>   fewer iterations than the cold solve
   subroutine check_uv3d_reload(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(model_moz_3d_parameters_type) :: params
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(model_moz_3d_type), target :: cold, warm, moved
      real(wp) :: e_cold, e_warm, e_moved

      call new_context(ctx, nthreads=0, verbosity=0)
      params%tolerance = 1.0e-8_wp
      call vv_settings(params%solvent)
      call new_memory_store(store, log)
      call new_cartesian_point_grid(grid, err, 32, 32, 32, 0.6_wp, margin=8.0_wp)
      if (.not. allocated(err)) call new_water_3d(ctx, grid, params, store, 0.0_wp, cold, err)
      e_cold = 0.0_wp
      if (.not. allocated(err)) call model_energy(cold, e_cold, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the cold solve did not save exactly its final field")
      if (allocated(error)) return
      associate (saved => log%saved(1))
         call check(error, saved%converged .and. saved%iterations > 0 .and. saved%kind == checkpoint_kind_uv3d &
            & .and. all(saved%grid%shape == [32, 32, 32]) .and. allocated(saved%grid%origin) &
            & .and. saved%grid%kind == "cartesian" .and. .not. allocated(saved%grid%xyz) &
            & .and. saved%tail == tail_known .and. saved%alpha > 0.0_wp .and. size(saved%solute_charges) == 3 &
            & .and. size(saved%solvent_charges) == 3, more="unexpected saved 3D UV state")
      end associate
      if (allocated(error)) return

      store%preset = log%saved(1)
      call new_water_3d(ctx, grid, params, store, 0.0_wp, warm, err)
      e_warm = 0.0_wp
      if (.not. allocated(err)) call model_energy(warm, e_warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the converged reload iterated")
      if (allocated(error)) return
      call check(error, e_warm, e_cold, thr=1.0e-12_wp, more="reloaded energy differs")
      if (allocated(error)) return

      call new_water_3d(ctx, grid, params, store, 0.05_wp, moved, err)
      e_moved = 0.0_wp
      if (.not. allocated(err)) call model_energy(moved, e_moved, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 2, more="the displaced solute did not iterate from the guess")
      if (allocated(error)) return
      call check(error, log%saved(2)%converged .and. log%saved(2)%iterations < log%saved(1)%iterations, &
         & more="the warm start was not faster than the cold solve")
   end subroutine check_uv3d_reload

   !> States on atom-centered grids seed only the same atomic grids on built-in Lebedev
   !> directions
   !>
   !> - Refused for a header expecting fields, with the reason: other shell radii, angular
   !>   points or angular directions, no point positions
   !> - Refused: decreasing positions, positions not matching the field, another column count
   !> - The same atomic grids give a guess
   !>
   !> @param[out] error  test failure
   subroutine check_points_compatible(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moz_checkpoint_state_type) :: state, header, other
      type(memory_store_type) :: store
      type(store_log_type), target :: log
      real(wp), allocatable :: guess(:, :)
      logical, allocatable :: matched(:)

      call new_context(ctx, nthreads=0, verbosity=0)
      call molecular_state(state, stored_points)
      call molecular_state(header, [2, 3, 4, 5, 8, 13, 21, 24])
      deallocate (header%tr)
      header%field_shape = [8, 2]
      call state%check_compatible(header, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call new_memory_store(store, log)
      store%header = header
      store%preset = state
      call store%load_guess(ctx, 1.0_wp, spread([0.0_wp, 0.0_wp], 1, 8), 4, guess, matched)
      call check(error, allocated(guess), more="the same atomic grids gave no guess")
      if (allocated(error)) return

      other = state
      other%grid%shell_radius(2) = 1.1_wp*other%grid%shell_radius(2)
      call offered_unused(other, "other shell radii")
      if (allocated(error)) return
      other = state
      other%grid%shell_points = [6, 6, 6, 14]
      call offered_unused(other, "other angular points")
      if (allocated(error)) return
      other = state
      deallocate (other%grid%angular_grid)
      call offered_unused(other, "a state without its angular grid")
      if (allocated(error)) return
      other = state
      other%grid%angular_grid = "custom"
      call offered_unused(other, "a custom angular grid")
      if (allocated(error)) return
      other = state
      deallocate (other%grid%point_index)
      call offered_unused(other, "a state without point positions")
      if (allocated(error)) return
      other = state
      other%grid%point_index(3) = 1
      call refused(other, "decreasing point positions", "point positions do not increase")
      if (allocated(error)) return
      other = state
      other%grid%point_index = other%grid%point_index(2:)
      call refused(other, "point positions not matching the field", "fields do not match the grid")
      if (allocated(error)) return
      other = state
      other%field_shape = [7, 3]
      call refused(other, "another column count", "field shape differs")
      if (allocated(error)) return
      call ctx%delete()
   contains
      !> Expect a refusal of a state against the header
      !>
      !> @param[in] candidate  state to check
      !> @param[in] label      what is wrong with it
      !> @param[in] expected   part of the expected refusal
      subroutine refused(candidate, label, expected)
         !> State to check
         type(moz_checkpoint_state_type), intent(in) :: candidate
         !> What is wrong, and part of the expected refusal
         character(len=*), intent(in) :: label, expected
         type(moist_error), allocatable :: err
         call candidate%check_compatible(header, err)
         call check(error, allocated(err), more=label//" accepted")
         if (allocated(error)) return
         call check(error, index(err%message, expected) > 0, more=label//": "//err%message)
      end subroutine refused

      !> Expect a state to be refused for the header, for want of a method, and to give no guess
      !>
      !> Other atomic grids go to the transfer by coordinates, which needs the solute
      !> coordinates these synthetic states do not have
      !>
      !> @param[in] candidate  state to check
      !> @param[in] label      how it differs
      subroutine offered_unused(candidate, label)
         !> State to check
         type(moz_checkpoint_state_type), intent(in) :: candidate
         !> How it differs
         character(len=*), intent(in) :: label
         type(moist_error), allocatable :: err
         call candidate%check_compatible(header, err)
         call check(error, allocated(err), more=label//" offered")
         if (allocated(error)) return
         call check(error, index(err%message, "solute coordinates") > 0, more=label//": "//err%message)
         if (allocated(error)) return
         store%preset = candidate
         call store%load_guess(ctx, 1.0_wp, spread([0.0_wp, 0.0_wp], 1, 8), 4, guess, matched)
         call check(error, .not. allocated(guess), more=label//" gave a guess")
      end subroutine offered_unused
   end subroutine check_points_compatible

   !> Molecular grids report built-in Lebedev points, but not a custom generator with the same counts
   !>
   !> @param[out] error  test failure
   subroutine check_points_angular(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: builtin, custom
      type(permuted_angular_type) :: angular
      type(structure_type) :: mol

      call new(mol, [8, 1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.43_wp, 1.1_wp, 0.0_wp, &
         & -1.43_wp, 1.1_wp, 0.0_wp], [3, 3]))
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(builtin, err, recipe=recipe, reciprocal=.false.)
      if (.not. allocated(err)) call builtin%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, builtin%has_lebedev_points(), more="built-in Lebedev points not recognized")
      if (allocated(error)) return
      deallocate (recipe%angular)
      allocate (recipe%angular, source=angular)
      call new_molecular_point_grid(custom, err, recipe=recipe, reciprocal=.false.)
      if (.not. allocated(err)) call custom%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. custom%has_lebedev_points(), more="custom angular directions taken as built-in")
   end subroutine check_points_angular

   !> A stored state seeds the points with the same position and marks them
   !>
   !> - Covered points get t - beta u_lr of their stored point
   !> - Without a mask, only a state covering every point is used
   !> - A state covering no point gives no guess
   !>
   !> @param[out] error  test failure
   subroutine check_points_guess(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(memory_store_type) :: store
      type(store_log_type), target :: log
      real(wp), parameter :: beta = 1.7_wp
      integer, parameter :: target_points(7) = [2, 3, 5, 8, 13, 21, 24]
      logical, parameter :: expected(7) = [.true., .false., .true., .true., .true., .false., .false.]
      real(wp), allocatable :: ur_lr(:, :), guess(:, :)
      logical, allocatable :: matched(:)
      integer :: i, j

      call new_context(ctx, nthreads=1, verbosity=0)
      call new_memory_store(store, log)
      call molecular_state(store%preset, stored_points)
      call molecular_state(store%header, target_points)
      deallocate (store%header%tr)
      allocate (ur_lr(size(target_points), 2))
      do j = 1, 2
         do i = 1, size(target_points)
            ur_lr(i, j) = 0.01_wp*real(i + 3*j, wp)
         end do
      end do

      call store%load_guess(ctx, beta, ur_lr, 4, guess, matched)
      call check(error, allocated(guess) .and. allocated(matched), more="partial state gave no guess")
      if (allocated(error)) return
      call check(error, all(matched .eqv. expected), more="wrong covered points")
      if (allocated(error)) return
      do i = 1, size(target_points)
         if (.not. expected(i)) cycle
         j = findloc(stored_points, target_points(i), dim=1)
         call check(error, all(guess(i, :) == store%preset%tr(j, :) - beta*ur_lr(i, :)), &
            & more="covered point did not get its stored value")
         if (allocated(error)) return
      end do

      call store%load_guess(ctx, beta, ur_lr, 4, guess)
      call check(error, .not. allocated(guess), more="partial state used without a mask")
      if (allocated(error)) return

      store%header%grid%point_index = [3, 4, 6, 9, 10, 11, 12]
      call store%load_guess(ctx, beta, ur_lr, 4, guess, matched)
      call check(error, .not. allocated(guess) .and. .not. allocated(matched), more="state covering no point used")
      if (allocated(error)) return

      store%header%grid%point_index = stored_points
      call store%load_guess(ctx, beta, ur_lr, 4, guess)
      call check(error, allocated(guess), more="state covering every point refused without a mask")
      call ctx%delete()
   end subroutine check_points_guess

   !> Point positions and atomic grids survive a file round trip; a file without point
   !> positions is still offered, with its coordinates
   !>
   !> @param[out] error  test failure
   subroutine check_file_points(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err, warning
      type(moz_file_checkpoint_type) :: store
      type(moz_checkpoint_state_type) :: state, loaded
      character(len=:), allocatable :: source
      logical :: found

      call molecular_state(state, stored_points)
      call new_file_checkpoint(store, file_points)
      call store%save(state, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      store%header = header_of(state)
      store%header%grid%point_index = [2, 3, 5]
      store%header%grid%shell_offset = state%grid%shell_offset
      store%header%grid%shell_radius = state%grid%shell_radius
      store%header%grid%shell_points = state%grid%shell_points
      call store%load(loaded, found, source, warning)
      call check(error, found .and. .not. allocated(warning), more="saved molecular state not found")
      if (allocated(error)) then
         call delete_file(file_points, cleanup_error)
         if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
         if (allocated(error)) return
         return
      end if
      call check(error, all(loaded%grid%point_index == state%grid%point_index) &
         & .and. all(loaded%grid%shell_offset == state%grid%shell_offset) &
         & .and. all(loaded%grid%shell_radius == state%grid%shell_radius) &
         & .and. all(loaded%grid%shell_points == state%grid%shell_points) &
         & .and. loaded%grid%angular_grid == state%grid%angular_grid .and. loaded%grid%kind == "points" &
         & .and. all(loaded%grid%xyz == state%grid%xyz) .and. all(loaded%grid%w == state%grid%w) &
         & .and. all(loaded%grid%owner == state%grid%owner) .and. all(loaded%grid%kmodes == state%grid%kmodes) &
         & .and. all(loaded%grid%kspacing == state%grid%kspacing) .and. all(loaded%grid%kref == state%grid%kref) &
         & .and. all(loaded%tr == state%tr), more="point positions or atomic grids changed")
      if (allocated(error)) then
         call delete_file(file_points, cleanup_error)
         if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
         if (allocated(error)) return
         return
      end if

      deallocate (state%grid%point_index)
      call store%save(state, err)
      if (.not. allocated(err)) call store%load(loaded, found, source, warning)
      call delete_file(file_points, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, found .and. .not. allocated(loaded%grid%point_index) .and. allocated(loaded%grid%xyz), &
         & more="state without point positions was not offered with its coordinates")
   end subroutine check_file_points

   !> A 3D UV model on a molecular grid describes its points in the saved states
   !>
   !> One iteration per solve; the second structure moves an atom, which changes the
   !> retained points, and starts from the first state
   !>
   !> @param[out] error  test failure
   subroutine check_uv3d_molecular(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(model_moz_3d_parameters_type) :: params
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: grid
      type(model_moz_3d_type), target :: first, second
      real(wp) :: energy
      logical :: same

      call new_context(ctx, nthreads=0, verbosity=0)
      params%maxiter = 1
      ! KH is the first rung of the closure ramp, so the single iteration is saved
      params%closure = "KH"
      call vv_settings(params%solvent)
      call new_memory_store(store, log)
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(grid, err, recipe=recipe, pruning_threshold=0.1_wp, &
         & dr=1.1_wp)
      if (.not. allocated(err)) call new_water_3d(ctx, grid, params, store, 0.0_wp, first, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      energy = 0.0_wp
      call model_energy(first, energy, err)
      call check(error, unconverged(err) .and. size(log%saved) == 1, more="the first solve did not save its field")
      if (allocated(error)) return
      select type (realized => first%grid)
      type is (moist_math_grid_3d_molecular_type)
         associate (saved => log%saved(1))
            call check(error, all(saved%grid%shape == [sum(realized%shell_npts)]) .and. saved%grid%spacing == 0.0_wp &
               & .and. all(saved%grid%point_index == realized%raw_index) .and. size(saved%tr, 1) == realized%ngrid &
               & .and. all(saved%grid%shell_offset == realized%atom_shell_offset) &
               & .and. all(saved%grid%shell_radius == realized%shell_r) &
               & .and. all(saved%grid%shell_points == realized%shell_npts) .and. saved%grid%angular_grid == "lebedev" &
               & .and. saved%grid%kind == "points" .and. all(saved%grid%xyz == realized%xyz) &
               & .and. all(saved%grid%w == realized%w) .and. all(saved%grid%owner == realized%owner) &
               & .and. (allocated(saved%grid%kmodes) .eqv. realized%has_reciprocal()), &
               & more="molecular grid not described")
         end associate
      class default
         call check(error, .false., more="the model grid is not molecular")
      end select
      if (allocated(error)) return

      store%preset = log%saved(1)
      call new_water_3d(ctx, grid, params, store, 0.3_wp, second, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call model_energy(second, energy, err)
      call check(error, unconverged(err) .and. size(log%saved) == 2, more="the second solve did not save its field")
      if (allocated(error)) return
      associate (a => log%saved(1)%grid%point_index, b => log%saved(2)%grid%point_index)
         same = size(a) == size(b)
         if (same) same = all(a == b)
         call check(error, .not. same .and. size(b) == size(log%saved(2)%tr, 1), &
            & more="the moved atom did not change the retained points")
      end associate
      call ctx%delete()
   contains
      !> Whether a solve failed only by running out of iterations
      !>
      !> @param[in] err  energy failure
      logical function unconverged(err)
         !> Energy failure
         type(moist_error), allocatable, intent(in) :: err
         unconverged = .false.
         if (allocated(err)) unconverged = index(err%message, "did not converge") > 0
      end function unconverged
   end subroutine check_uv3d_molecular

   !> A memory store saves into the slot of the caller and offers it back
   !>
   !> - Copies of the store share the slot, and each save replaces it
   !> - An empty slot offers nothing; a slot of another grid is offered for transfer
   !> - A slot filled by hand needs no field shape
   !> - A store without a slot saves and offers nothing
   !>
   !> @param[out] error  test failure
   subroutine check_memory_store(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moz_checkpoint_state_type), target :: slot
      type(moz_checkpoint_state_type) :: state, loaded
      type(moz_memory_checkpoint_type) :: store, unset
      class(moz_checkpoint_type), allocatable :: copy
      character(len=:), allocatable :: source
      logical :: found

      call synthetic_state(state)
      call new_memory_checkpoint(store, slot)
      store%header = header_of(state)
      call check(error, store%describe() == "memory", more="memory store not described")
      if (allocated(error)) return
      call store%load(loaded, found, source, err)
      call check(error, .not. (found .or. allocated(err)), more="an empty slot offered a state")
      if (allocated(error)) return

      allocate (copy, source=store)
      call copy%save(state, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(slot%tr), more="a copy of the store did not save into the slot")
      if (allocated(error)) return
      state%tr = 2.0_wp*state%tr
      call copy%save(state, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(slot%tr == state%tr), more="a second save did not replace the slot")
      if (allocated(error)) return
      call store%load(loaded, found, source, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, found .and. source == "memory", more="the saved slot was not offered")
      if (allocated(error)) return
      call check(error, all(loaded%tr == state%tr), more="the offered state differs from the saved one")
      if (allocated(error)) return

      slot = state
      if (allocated(slot%field_shape)) deallocate (slot%field_shape)
      call store%load(loaded, found, source, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, found, more="a slot without field shape was not offered")
      if (allocated(error)) return
      store%header%grid%shape = [2, 2, 4]
      call store%load(loaded, found, source, err)
      call check(error, found .and. .not. allocated(err), more="a slot of another grid was not offered")
      if (allocated(error)) return
      store%header%solute_sites = [9]
      call store%load(loaded, found, source, err)
      call check(error, .not. found .and. allocated(err), more="a slot of another solute was offered")
      if (allocated(error)) return

      call unset%save(state, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      unset%header = header_of(state)
      call unset%load(loaded, found, source, err)
      call check(error, .not. (found .or. allocated(err)) .and. len(unset%describe()) == 0, &
         & more="a store without a slot offered a state")
   end subroutine check_memory_store

   !> A memory store carries a 3D UV solution across a geometry update of one model
   !>
   !> - The cold solve leaves its converged field in the slot of the caller
   !> - After the solute moves, the same model starts from the slot, converges in
   !>   fewer iterations and replaces the slot with its own field
   !> - A new model at the moved structure accepts the slot without iterating
   !>
   !> @param[out] error  test failure
   subroutine check_memory_update(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moz_checkpoint_state_type), target :: slot
      type(moz_memory_checkpoint_type) :: store
      type(model_moz_3d_parameters_type) :: params
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(model_moz_3d_type), target :: model, fresh
      type(solvent_vv_type) :: solvent
      type(structure_type) :: mol
      real(wp) :: energy, e_fresh
      integer :: cold_iterations

      call new_context(ctx, nthreads=0, verbosity=0)
      params%tolerance = 1.0e-8_wp
      call vv_settings(params%solvent)
      call new_memory_checkpoint(store, slot)
      call new_cartesian_point_grid(grid, err, 32, 32, 32, 0.6_wp, margin=8.0_wp)
      if (.not. allocated(err)) call new_water_3d(ctx, grid, params, store, 0.0_wp, model, err)
      energy = 0.0_wp
      if (.not. allocated(err)) call model_energy(model, energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(slot%tr), more="the cold solve left no field in the slot")
      if (allocated(error)) return
      call check(error, slot%converged .and. slot%iterations > 0 .and. slot%kind == checkpoint_kind_uv3d, &
         & more="unexpected state in the slot")
      if (allocated(error)) return
      cold_iterations = slot%iterations

      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) then
         mol = solvent%structure()
         mol%xyz(1, 1) = mol%xyz(1, 1) + 0.05_wp
         call model%update(mol, err)
      end if
      energy = 0.0_wp
      if (.not. allocated(err)) call model_energy(model, energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, slot%converged .and. slot%iterations > 0 .and. slot%iterations < cold_iterations, &
         & more="the moved solute did not start from the slot")
      if (allocated(error)) return
      call check(error, slot%solute_xyz(1, 1), mol%xyz(1, 1), thr=1.0e-12_wp, &
         & more="the slot does not hold the solve of the moved solute")
      if (allocated(error)) return

      ! Provenance only; a save would overwrite it
      slot%iterations = -1
      call new_water_3d(ctx, grid, params, store, 0.05_wp, fresh, err)
      e_fresh = 0.0_wp
      if (.not. allocated(err)) call model_energy(fresh, e_fresh, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, slot%iterations == -1, more="the converged slot was solved again")
      if (allocated(error)) return
      call check(error, e_fresh, energy, thr=1.0e-12_wp, more="energy from the slot differs")
   end subroutine check_memory_update

   !> A memory store reloads a converged VV solution without iterations
   !>
   !> - The cold solve leaves its converged field at full coupling in the slot
   !> - A new solvent accepts the slot without iterating and has the same h(r) and mu
   !>
   !> @param[out] error  test failure
   subroutine check_memory_vv(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moz_checkpoint_state_type), target :: slot
      type(moz_memory_checkpoint_type) :: store
      type(solvent_vv_parameters_type) :: params
      type(solvent_vv_type) :: cold, warm
      real(wp) :: mu_cold, mu_warm, dh

      call new_context(ctx, nthreads=0, verbosity=0)
      call vv_settings(params)
      call new_memory_checkpoint(store, slot)
      call solve_vv(ctx, params, store, cold, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(slot%tr), more="the cold solve left no field in the slot")
      if (allocated(error)) return
      call check(error, slot%converged .and. slot%iterations > 0 .and. slot%lambda == 1.0_wp &
         & .and. slot%kind == checkpoint_kind_vv, more="unexpected state in the slot")
      if (allocated(error)) return

      ! Provenance only; a save would overwrite it
      slot%iterations = -1
      call solve_vv(ctx, params, store, warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, slot%iterations == -1 .and. warm%is_solved(params), &
         & more="the converged slot was solved again or left the solvent unsolved")
      if (allocated(error)) return
      call vv_compare(cold, warm, mu_cold, mu_warm, dh, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, dh < 1.0e-12_wp, more="h(r) from the slot differs")
      if (allocated(error)) return
      call check(error, mu_warm, mu_cold, thr=1.0e-13_wp, more="chemical potential from the slot differs")
   end subroutine check_memory_vv

   !> A memory store carries a 1D UV solution across a geometry update of one model
   !>
   !> - The cold solve leaves its converged field in the slot
   !> - After the solute moves, the same model starts from the slot, converges in
   !>   fewer iterations and replaces the slot with its own field
   !> - Another model updated to the moved structure accepts the slot without iterating
   !>
   !> @param[out] error  test failure
   subroutine check_memory_uv1d(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(moz_checkpoint_state_type), target :: slot
      type(moz_memory_checkpoint_type) :: store
      type(model_moz_1d_parameters_type) :: params
      type(model_moz_1d_type), target :: model, fresh
      type(solvent_vv_type) :: solvent
      type(structure_type) :: mol
      real(wp) :: energy, e_fresh
      integer :: cold_iterations

      call new_context(ctx, nthreads=0, verbosity=0)
      params%tolerance = 1.0e-8_wp
      call vv_settings(params%solvent)
      call new_memory_checkpoint(store, slot)
      call new_water_1d(ctx, params, store, model, err)
      energy = 0.0_wp
      if (.not. allocated(err)) call model_energy(model, energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(slot%tr), more="the cold solve left no field in the slot")
      if (allocated(error)) return
      call check(error, slot%converged .and. slot%iterations > 0 .and. slot%kind == checkpoint_kind_uv1d, &
         & more="unexpected state in the slot")
      if (allocated(error)) return
      cold_iterations = slot%iterations

      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) then
         mol = solvent%structure()
         mol%xyz(1, 1) = mol%xyz(1, 1) + 0.05_wp
         call model%update(mol, err)
      end if
      energy = 0.0_wp
      if (.not. allocated(err)) call model_energy(model, energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, slot%converged .and. slot%iterations > 0 .and. slot%iterations < cold_iterations, &
         & more="the moved solute did not start from the slot")
      if (allocated(error)) return
      call check(error, slot%solute_xyz(1, 1), mol%xyz(1, 1), thr=1.0e-12_wp, &
         & more="the slot does not hold the solve of the moved solute")
      if (allocated(error)) return

      ! Provenance only; a save would overwrite it
      slot%iterations = -1
      call new_water_1d(ctx, params, store, fresh, err)
      if (.not. allocated(err)) call fresh%update(mol, err)
      e_fresh = 0.0_wp
      if (.not. allocated(err)) call model_energy(fresh, e_fresh, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, slot%iterations == -1, more="the converged slot was solved again")
      if (allocated(error)) return
      call check(error, e_fresh, energy, thr=1.0e-12_wp, more="energy from the slot differs")
   end subroutine check_memory_uv1d

   !> Small grids and tolerance of the VV solves
   !>
   !> @param[out] params  settings
   subroutine vv_settings(params)
      !> Settings
      type(solvent_vv_parameters_type), intent(out) :: params
      params%closure = "KH"
      params%npts = 1024
      params%dr = 0.06_wp
      params%tolerance = 1.0e-9_wp
   end subroutine vv_settings

   !> Empty memory store with its log
   !>
   !> @param[out]   store  store
   !> @param[inout] log    log receiving the saves
   subroutine new_memory_store(store, log)
      !> Store
      type(memory_store_type), intent(out) :: store
      !> Log of the saves
      type(store_log_type), intent(inout), target :: log
      if (allocated(log%saved)) deallocate (log%saved)
      allocate (log%saved(0))
      store%log => log
   end subroutine new_memory_store

   !> Solve cSPC/E water with a store
   !>
   !> @param[in]  ctx      run context
   !> @param[in]  params   settings
   !> @param[in]  store    store, copied into the solvent
   !> @param[out] solvent  solved solvent
   !> @param[out] err      setup or solve failure
   subroutine solve_vv(ctx, params, store, solvent, err)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Settings
      type(solvent_vv_parameters_type), intent(in) :: params
      !> Store
      class(moz_checkpoint_type), intent(in) :: store
      !> Solved solvent
      type(solvent_vv_type), intent(out) :: solvent
      !> Setup or solve failure
      type(moist_error), allocatable, intent(out) :: err
      call new_cspce_water_solvent(ctx, solvent, err)
      if (allocated(err)) return
      call solvent%set_checkpoint(store)
      call solvent%solve(err, params)
   end subroutine solve_vv

   !> Chemical potentials and the largest h(r) difference of two solved solvents
   !>
   !> @param[in]  a     first solvent
   !> @param[in]  b     second solvent
   !> @param[out] mu_a  chemical potential of a, Hartree
   !> @param[out] mu_b  chemical potential of b, Hartree
   !> @param[out] dh    largest absolute h(r) difference
   !> @param[out] err   unsolved solvent
   subroutine vv_compare(a, b, mu_a, mu_b, dh, err)
      !> First solvent
      type(solvent_vv_type), intent(in) :: a
      !> Second solvent
      type(solvent_vv_type), intent(in) :: b
      !> Chemical potentials
      real(wp), intent(out) :: mu_a, mu_b
      !> Largest h(r) difference
      real(wp), intent(out) :: dh
      !> Unsolved solvent
      type(moist_error), allocatable, intent(out) :: err
      type(solvent_vv_thermodynamics_type) :: ta, tb
      real(wp), allocatable :: cr(:, :), tr(:, :), ha(:, :), hb(:, :)

      mu_a = 0.0_wp
      mu_b = 0.0_wp
      dh = huge(1.0_wp)
      call a%thermodynamics(ta, err)
      if (.not. allocated(err)) call b%thermodynamics(tb, err)
      if (.not. allocated(err)) call a%correlations(cr, tr, ha, err)
      if (.not. allocated(err)) call b%correlations(cr, tr, hb, err)
      if (allocated(err)) return
      mu_a = ta%chemical_potential
      mu_b = tb%chemical_potential
      dh = maxval(abs(ha - hb))
   end subroutine vv_compare

   !> Updated 1D model of a water solute in cSPC/E water with a store
   !>
   !> @param[in]  ctx     run context
   !> @param[in]  params  settings
   !> @param[in]  store   UV store, copied into the model
   !> @param[out] model   updated model
   !> @param[out] err     construction or update failure
   subroutine new_water_1d(ctx, params, store, model, err)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Settings
      type(model_moz_1d_parameters_type), intent(in) :: params
      !> UV store
      class(moz_checkpoint_type), intent(in) :: store
      !> Updated model
      type(model_moz_1d_type), intent(out) :: model
      !> Construction or update failure
      type(moist_error), allocatable, intent(out) :: err
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_radial_type) :: rgrid, kgrid

      call new_cspce_water_solvent(ctx, solvent, err)
      if (.not. allocated(err)) call new_uniform_radial_pair(rgrid, kgrid, params%solvent%npts, params%solvent%dr, err)
      if (.not. allocated(err)) call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, err, params)
      if (.not. allocated(err)) call model%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call model%potential%add(coulomb_spce, err)
      if (.not. allocated(err)) call model%update(solvent%structure(), err)
      if (.not. allocated(err)) call model%set_checkpoint(store)
   end subroutine new_water_1d

   !> Updated 3D model of a water solute in cSPC/E water with a store
   !>
   !> @param[in]  ctx     run context
   !> @param[in]  grid    grid template
   !> @param[in]  params  settings
   !> @param[in]  store   UV store, copied into the model
   !> @param[in]  shift   displacement of the first solute atom along x, bohr
   !> @param[out] model   updated model
   !> @param[out] err     construction or update failure
   subroutine new_water_3d(ctx, grid, params, store, shift, model, err)
      !> Run context
      type(moist_context_type), intent(in), target :: ctx
      !> Grid template
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Settings
      type(model_moz_3d_parameters_type), intent(in) :: params
      !> UV store
      class(moz_checkpoint_type), intent(in) :: store
      !> Displacement, bohr
      real(wp), intent(in) :: shift
      !> Updated model
      type(model_moz_3d_type), intent(out) :: model
      !> Construction or update failure
      type(moist_error), allocatable, intent(out) :: err
      type(solvent_vv_type) :: solvent
      type(structure_type) :: mol

      call new_cspce_water_solvent(ctx, solvent, err)
      if (allocated(err)) return
      mol = solvent%structure()
      mol%xyz(1, 1) = mol%xyz(1, 1) + shift
      call new_moz_3d_model(model, ctx, grid, solvent, err, params)
      if (.not. allocated(err)) call model%potential%add(lj_cspce, err)
      if (.not. allocated(err)) call model%potential%add(coulomb_spce, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (.not. allocated(err)) call model%set_checkpoint(store)
   end subroutine new_water_3d

   !> Energy through a coupling of the model
   !>
   !> @param[inout] model   updated model
   !> @param[inout] energy  energy, Hartree
   !> @param[out]   err     coupling or energy failure
   subroutine model_energy(model, energy, err)
      !> Model
      class(solvation_model_type), intent(inout), target :: model
      !> Energy
      real(wp), intent(inout) :: energy
      !> Coupling or energy failure
      type(moist_error), allocatable, intent(out) :: err
      type(coupling_type), pointer :: coupling
      call model%new_coupling(coupling, err)
      if (allocated(err)) return
      call model%prepare_energy(coupling, err)
      if (.not. allocated(err)) call model%get_energy(coupling, energy, err)
      call model%release_coupling(coupling)
   end subroutine model_energy

   !> A VV state of a coarser grid seeds a finer one through the radial spline
   !>
   !> - Every node of the finer grid but the last lies within the stored range
   !> - The warm start is kept (no rejection, no cold retry) and reaches the chemical
   !>   potential and h(r) of a cold solve
   !>
   !> @param[out] error  test failure
   subroutine check_transfer_vv(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log, quiet_log
      type(memory_store_type) :: store, quiet
      type(solvent_vv_parameters_type) :: coarse, fine
      type(solvent_vv_type) :: source, cold, warm
      real(wp) :: mu_cold, mu_warm, dh
      logical :: splined
      character(len=*), parameter :: logfile = "moist-test-checkpoint-transfer-vv.log"

      call vv_settings(coarse)
      fine = coarse
      fine%npts = 2*coarse%npts
      fine%dr = 0.5_wp*coarse%dr
      call new_memory_store(store, log)
      call new_memory_store(quiet, quiet_log)
      call new_context(ctx, nthreads=0, verbosity=0)
      call solve_vv(ctx, coarse, store, source, err)
      if (.not. allocated(err)) call solve_vv(ctx, fine, quiet, cold, err)
      call ctx%delete()
      call check_moist_error(error, err)
      if (allocated(error)) return
      store%preset = log%saved(size(log%saved))
      call new_context(ctx, nthreads=0, verbosity=2, logfile=logfile)
      call solve_vv(ctx, fine, store, warm, err)
      call ctx%delete()
      if (.not. allocated(err)) call vv_compare(cold, warm, mu_cold, mu_warm, dh, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      splined = logged(logfile, "spline, 2047 of 2048 points") .and. warm_start(logfile)
      call check(error, size(log%saved) == 2 .and. size(quiet_log%saved) == 1, &
         & more="expected exactly the source, warm and cold reference saves")
      if (allocated(error)) return
      call check(error, log%saved(2)%iterations < quiet_log%saved(1)%iterations, &
         & more="the transferred guess did not reduce the cold iteration count")
      if (allocated(error)) return
      call delete_file(logfile, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, splined, more="the coarser VV state was not splined or not kept")
      if (allocated(error)) return
      call check(error, mu_warm, mu_cold, rel=.true., thr=1.0e-7_wp, more="warm VV chemical potential differs")
      if (allocated(error)) return
      call check(error, dh < 1.0e-6_wp, more="warm VV h(r) differs")
   end subroutine check_transfer_vv

   !> A 1D UV state seeds a longer radial grid of the same spacing: the shared nodes are
   !> copied, nodes beyond the stored range start cold; the warm start is kept and reaches
   !> the energy of a cold solve
   !>
   !> @param[out] error  test failure
   subroutine check_transfer_uv1d(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log, quiet_log
      type(memory_store_type) :: store, quiet
      type(model_moz_1d_parameters_type) :: short, long
      type(model_moz_1d_type) :: source, cold, warm
      real(wp) :: e_cold, e_warm
      logical :: partial
      character(len=*), parameter :: logfile = "moist-test-checkpoint-transfer-uv1d.log"

      call vv_settings(short%solvent)
      short%closure = "KH"
      short%tolerance = 1.0e-9_wp
      long = short
      long%solvent%npts = 2*short%solvent%npts
      call new_memory_store(store, log)
      call new_memory_store(quiet, quiet_log)
      call new_context(ctx, nthreads=0, verbosity=0)
      e_cold = 0.0_wp
      call new_water_1d(ctx, short, store, source, err)
      if (.not. allocated(err)) call model_energy(source, e_cold, err)
      e_cold = 0.0_wp
      if (.not. allocated(err)) call new_water_1d(ctx, long, quiet, cold, err)
      if (.not. allocated(err)) call model_energy(cold, e_cold, err)
      call ctx%delete()
      call check_moist_error(error, err)
      if (allocated(error)) return
      store%preset = log%saved(size(log%saved))
      call new_context(ctx, nthreads=0, verbosity=4, logfile=logfile)
      e_warm = 0.0_wp
      call new_water_1d(ctx, long, store, warm, err)
      if (.not. allocated(err)) call model_energy(warm, e_warm, err)
      call ctx%delete()
      call check_moist_error(error, err)
      if (allocated(error)) return
      partial = logged(logfile, "copy, ") .and. logged(logfile, "the others cold") .and. warm_start(logfile)
      call check(error, size(log%saved) == 2 .and. size(quiet_log%saved) == 1, &
         & more="expected exactly the source, warm and cold reference saves")
      if (allocated(error)) return
      call check(error, log%saved(2)%iterations < quiet_log%saved(1)%iterations, &
         & more="the transferred guess did not reduce the cold iteration count")
      if (allocated(error)) return
      call delete_file(logfile, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, partial, more="the shorter 1D UV state was not copied onto part of the longer grid")
      if (allocated(error)) return
      call check(error, e_warm, e_cold, thr=1.0e-8_wp, more="warm 1D UV energy differs")
   end subroutine check_transfer_uv1d

   !> A 3D UV state seeds a larger box on the same lattice by copying the shared points
   !>
   !> @param[out] error  test failure
   subroutine check_transfer_uv3d_box(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: small, large
      type(moist_error), allocatable :: err

      call new_cartesian_point_grid(small, err, 24, 24, 24, 0.6_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_cartesian_point_grid(large, err, 32, 32, 32, 0.6_wp, margin=0.0_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call transfer_uv3d(error, small, large, 0.0_wp, "copy, ", "moist-test-checkpoint-transfer-box.log")
   end subroutine check_transfer_uv3d_box

   !> A 3D UV state seeds a finer lattice of the same box through its Fourier modes
   !>
   !> @param[out] error  test failure
   subroutine check_transfer_uv3d_finer(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: coarse, fine
      type(moist_error), allocatable :: err

      call new_cartesian_point_grid(coarse, err, 24, 24, 24, 0.6_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_cartesian_point_grid(fine, err, 36, 36, 36, 0.4_wp, margin=0.0_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call transfer_uv3d(error, coarse, fine, 0.0_wp, "fourier)", "moist-test-checkpoint-transfer-finer.log")
   end subroutine check_transfer_uv3d_finer

   !> A 3D UV state seeds the same lattice around a rotated solute through the rigid fit
   !>
   !> @param[out] error  test failure
   subroutine check_transfer_uv3d_rotated(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(moist_error), allocatable :: err

      call new_cartesian_point_grid(grid, err, 28, 28, 28, 0.6_wp, margin=0.0_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call transfer_uv3d(error, grid, grid, 0.6_wp, "fourier", "moist-test-checkpoint-transfer-rotated.log")
   end subroutine check_transfer_uv3d_rotated

   !> A Cartesian checkpoint seeds a molecular-grid solve through trilinear transfer
   !>
   !> @param[out] error  test failure
   subroutine check_transfer_cartesian_molecular(error)
      type(error_type), allocatable, intent(out) :: error
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(model_moz_3d_parameters_type) :: params
      type(model_moz_3d_type), target :: source, target
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(moist_math_grid_3d_molecular_type) :: molecular
      type(moist_math_grid_atomic_recipe_type) :: recipe
      real(wp) :: energy
      logical :: transferred
      character(len=*), parameter :: logfile = "moist-test-checkpoint-cartesian-molecular.log"

      call new_context(ctx, nthreads=0, verbosity=0)
      call vv_settings(params%solvent)
      params%closure = "KH"
      call new_memory_store(store, log)
      call new_cartesian_point_grid(cart, err, 24, 24, 24, 0.5_wp, margin=0.0_wp)
      if (.not. allocated(err)) call new_water_3d(ctx, cart, params, store, 0.0_wp, source, err)
      energy = 0.0_wp
      if (.not. allocated(err)) call model_energy(source, energy, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the Cartesian reference checkpoint was not saved")
      if (allocated(error)) return
      store%preset = log%saved(1)
      call ctx%delete()
      call new_context(ctx, nthreads=0, verbosity=4, logfile=logfile)
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(molecular, err, recipe=recipe, pruning_threshold=0.1_wp, &
         & dr=1.1_wp)
      ! This intentionally coarse molecular grid exercises seeding for one step.
      params%maxiter = 1
      if (.not. allocated(err)) call new_water_3d(ctx, molecular, params, store, 0.0_wp, target, err)
      if (.not. allocated(err)) call model_energy(target, energy, err)
      call ctx%delete()
      call check(error, allocated(err), more="the one-step molecular solve unexpectedly converged")
      if (allocated(error)) return
      call check(error, index(err%message, "3D UV solve did not converge") > 0, more=err%message)
      if (allocated(error)) return
      transferred = logged(logfile, "initial guess: memory") .and. logged(logfile, "trilinear") .and. &
         & logged(logfile, "Initial guess: RMS residual") .and. .not. logged(logfile, "ignoring guess")
      call delete_file(logfile, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, transferred, more="the Cartesian checkpoint was not accepted on the molecular grid")
   end subroutine check_transfer_cartesian_molecular

   !> A 3D UV state on a molecular grid seeds a Cartesian grid by moving least squares
   !>
   !> @param[out] error  test failure
   subroutine check_transfer_molecular(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_3d_cartesian_type) :: cart

      call new_cartesian_point_grid(cart, err, 16, 16, 16, 0.8_wp, margin=0.0_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call seed_from_molecular(error, cart, "moist-test-checkpoint-transfer-molecular.log")
   end subroutine check_transfer_molecular

   !> A 3D UV state on a molecular grid seeds a molecular grid of other radial shells by
   !> moving least squares
   !>
   !> @param[out] error  test failure
   subroutine check_transfer_uv3d_points(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: other

      call get_qc_handymod_recipe(recipe, err, nrad=10, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(other, err, recipe=recipe, pruning_threshold=0.1_wp, &
         & dr=1.1_wp)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call seed_from_molecular(error, other, "moist-test-checkpoint-transfer-points.log")
   end subroutine check_transfer_uv3d_points

   !> Save one iteration of water in water on a small molecular grid, then seed a solve on
   !> another grid: the state is moved by moving least squares and its residual is accepted
   !>
   !> Small molecular grids do not converge, so both solves stop after one iteration and
   !> the energies are not compared
   !>
   !> @param[out] error    test failure
   !> @param[in]  grid     grid of the seeded solve
   !> @param[in]  logfile  output of the seeded solve
   subroutine seed_from_molecular(error, grid, logfile)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Grid of the seeded solve
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Output of the seeded solve
      character(len=*), intent(in) :: logfile
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log
      type(memory_store_type) :: store
      type(model_moz_3d_parameters_type) :: params
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_3d_molecular_type) :: molecular
      type(model_moz_3d_type), target :: source, target
      real(wp) :: energy
      logical :: moved

      params%closure = "KH"
      params%maxiter = 1
      call vv_settings(params%solvent)
      call new_memory_store(store, log)
      call new_context(ctx, nthreads=0, verbosity=0)
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=5.0_wp)
      if (.not. allocated(err)) call new_molecular_point_grid(molecular, err, recipe=recipe, pruning_threshold=0.1_wp, &
         & dr=1.1_wp)
      if (.not. allocated(err)) call new_water_3d(ctx, molecular, params, store, 0.0_wp, source, err)
      energy = 0.0_wp
      if (.not. allocated(err)) call model_energy(source, energy, err)
      call ctx%delete()
      call check(error, allocated(err), more="the one-step 3D UV solve unexpectedly converged")
      if (allocated(error)) return
      call check(error, index(err%message, "3D UV solve did not converge") > 0, more=err%message)
      if (allocated(error)) return
      call check(error, size(log%saved) == 1, more="the molecular solve saved nothing")
      if (allocated(error)) return
      store%preset = log%saved(1)
      call new_context(ctx, nthreads=0, verbosity=4, logfile=logfile)
      call new_water_3d(ctx, grid, params, store, 0.0_wp, target, err)
      if (.not. allocated(err)) call model_energy(target, energy, err)
      call ctx%delete()
      call check(error, allocated(err), more="the one-step 3D UV solve unexpectedly converged")
      if (allocated(error)) return
      call check(error, index(err%message, "3D UV solve did not converge") > 0, more=err%message)
      if (allocated(error)) return
      moved = logged(logfile, "initial guess: memory") .and. logged(logfile, ", points") .and. &
         & .not. logged(logfile, "ignoring guess")
      call delete_file(logfile, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, moved, more="the molecular state was not moved, or its guess was refused")
   end subroutine seed_from_molecular

   !> Solve water in water on one grid, reload the state on another, compare with a cold solve
   !>
   !> - The source message names the method; the warm start is kept (no rejection, no cold retry)
   !>
   !> @param[out] error    test failure
   !> @param[in]  first    grid of the stored solve
   !> @param[in]  second   grid of the reloaded solve
   !> @param[in]  angle    rotation of the second solute about z, rad
   !> @param[in]  method   text the source message must hold
   !> @param[in]  logfile  output of the reloaded solve
   subroutine transfer_uv3d(error, first, second, angle, method, logfile)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Grid of the stored solve
      class(moist_math_grid_3d_type), intent(in) :: first
      !> Grid of the reloaded solve
      class(moist_math_grid_3d_type), intent(in) :: second
      !> Rotation of the second solute
      real(wp), intent(in) :: angle
      !> Expected method text
      character(len=*), intent(in) :: method
      !> Output of the reloaded solve
      character(len=*), intent(in) :: logfile
      !> File-cleanup failure, separate from the result under test
      type(moist_error), allocatable :: cleanup_error
      type(moist_error), allocatable :: err
      type(moist_context_type), target :: ctx
      type(store_log_type), target :: log, quiet_log
      type(memory_store_type) :: store, quiet
      type(model_moz_3d_parameters_type) :: params
      type(model_moz_3d_type), target :: source, cold, warm
      real(wp) :: e_cold, e_warm, rot(3, 3)
      type(structure_type) :: mol
      type(solvent_vv_type) :: water
      logical :: moved

      params%closure = "KH"
      params%tolerance = 1.0e-8_wp
      call vv_settings(params%solvent)
      call new_memory_store(store, log)
      call new_memory_store(quiet, quiet_log)
      call new_context(ctx, nthreads=0, verbosity=0)
      e_cold = 0.0_wp
      call new_water_3d(ctx, first, params, store, 0.0_wp, source, err)
      if (.not. allocated(err)) call model_energy(source, e_cold, err)
      e_cold = 0.0_wp
      if (.not. allocated(err)) call new_water_3d(ctx, second, params, quiet, 0.0_wp, cold, err)
      rot = reshape([cos(angle), sin(angle), 0.0_wp, -sin(angle), cos(angle), 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      if (.not. allocated(err)) call new_cspce_water_solvent(ctx, water, err)
      if (.not. allocated(err)) then
         mol = water%structure()
         mol%xyz = matmul(rot, mol%xyz)
         if (angle /= 0.0_wp) call cold%update(mol, err)
      end if
      if (.not. allocated(err)) call model_energy(cold, e_cold, err)
      call ctx%delete()
      call check_moist_error(error, err)
      if (allocated(error)) return
      store%preset = log%saved(size(log%saved))
      call new_context(ctx, nthreads=0, verbosity=4, logfile=logfile)
      e_warm = 0.0_wp
      call new_water_3d(ctx, second, params, store, 0.0_wp, warm, err)
      if (.not. allocated(err) .and. angle /= 0.0_wp) call warm%update(mol, err)
      if (.not. allocated(err)) call model_energy(warm, e_warm, err)
      call ctx%delete()
      call check_moist_error(error, err)
      if (allocated(error)) return
      moved = logged(logfile, "initial guess: memory") .and. logged(logfile, method) .and. warm_start(logfile)
      call check(error, size(log%saved) == 2 .and. size(quiet_log%saved) == 1, &
         & more="expected exactly the source, warm and cold reference saves")
      if (allocated(error)) return
      call check(error, log%saved(2)%iterations < quiet_log%saved(1)%iterations, &
         & more="the transferred guess did not reduce the cold iteration count")
      if (allocated(error)) return
      call delete_file(logfile, cleanup_error)
      if (allocated(cleanup_error)) call check_moist_error(error, cleanup_error)
      if (allocated(error)) return
      call check(error, moved, more="the stored state was not moved by "//method//" or not kept")
      if (allocated(error)) return
      call check(error, e_warm, e_cold, thr=1.0e-7_wp, more="warm 3D UV energy differs")
   end subroutine transfer_uv3d

   !> Whether the warm start of a solve was kept: no rejected guess and no cold retry
   !>
   !> @param[in] file  output file
   logical function warm_start(file)
      !> Output file
      character(len=*), intent(in) :: file
      warm_start = .not. (logged(file, "ignoring guess") .or. logged(file, "Retry"))
   end function warm_start

   !> Whether an output file holds a text
   !>
   !> @param[in] file  output file
   !> @param[in] text  text to look for
   logical function logged(file, text)
      !> Output file
      character(len=*), intent(in) :: file
      !> Text to look for
      character(len=*), intent(in) :: text
      character(len=512) :: line
      integer :: unit, stat

      logged = .false.
      open (newunit=unit, file=file, status="old", action="read", iostat=stat)
      if (stat /= 0) return
      do
         read (unit, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (index(line, text) > 0) logged = .true.
      end do
      close (unit)
   end function logged

   !> Synthetic 3D UV state with every optional dataset
   !>
   !> @param[out] state  state with fields (8, 2)
   subroutine synthetic_state(state)
      !> State
      type(moz_checkpoint_state_type), intent(out) :: state
      integer :: i, j

      state%kind = checkpoint_kind_uv3d
      state%grid%kind = "cartesian"
      state%grid%shape = [2, 2, 2]
      state%grid%spacing = 0.5_wp
      state%grid%origin = [-1.0_wp, -2.0_wp, -3.0_wp]
      state%solvent_sites = [8, 1]
      state%solute_sites = [17]
      state%solute_xyz = reshape([0.1_wp, 0.2_wp, 0.3_wp], [3, 1])
      state%tail = tail_known
      state%alpha = 0.8_wp
      state%solute_charges = [-1.0_wp]
      state%solvent_charges = [-0.8476_wp, 0.4238_wp]
      state%closure = "KH"
      state%theory = "RISM"
      state%temperature = 298.15_wp
      state%residual = 3.5e-9_wp
      state%lambda = 1.0_wp
      state%iterations = 17
      state%converged = .true.
      allocate (state%tr(8, 2))
      do j = 1, 2
         do i = 1, 8
            state%tr(i, j) = 0.1_wp*real(i, wp) - real(j, wp)
         end do
      end do
   end subroutine synthetic_state

   !> Synthetic 3D UV state on an atom-centered grid with 24 unpruned points
   !>
   !> Two atoms with two shells of six points each
   !>
   !> @param[out] state   state with fields (npoint, 2)
   !> @param[in]  points  increasing point positions
   subroutine molecular_state(state, points)
      !> State
      type(moz_checkpoint_state_type), intent(out) :: state
      !> Point positions
      integer, intent(in) :: points(:)
      integer :: i

      state%kind = checkpoint_kind_uv3d
      state%grid%kind = "points"
      state%grid%shape = [24]
      state%grid%spacing = 0.0_wp
      state%solvent_sites = [8, 1]
      state%solute_sites = [8, 1]
      state%grid%point_index = points
      allocate (state%grid%xyz(3, size(points)), state%grid%w(size(points)), state%grid%owner(size(points)))
      do i = 1, size(points)
         state%grid%xyz(:, i) = [0.1_wp*real(points(i), wp), -0.2_wp, 0.05_wp*real(i, wp)]
         state%grid%w(i) = 0.01_wp*real(i, wp)
         state%grid%owner(i) = merge(1, 2, points(i) <= 12)
      end do
      state%grid%kmodes = [4, 6, 8]
      state%grid%kspacing = [0.5_wp, 0.4_wp, 0.3_wp]
      state%grid%kref = [0.1_wp, 0.0_wp, -0.1_wp]
      state%grid%shell_offset = [1, 3, 5]
      state%grid%shell_radius = [0.5_wp, 1.5_wp, 0.4_wp, 1.2_wp]
      state%grid%shell_points = [6, 6, 6, 6]
      state%grid%angular_grid = "lebedev"
      state%closure = "KH"
      state%temperature = 298.15_wp
      state%field_shape = [size(points), 2]
      allocate (state%tr(size(points), 2))
      do i = 1, size(points)
         state%tr(i, :) = [real(points(i), wp), -0.5_wp*real(points(i), wp)]
      end do
   end subroutine molecular_state

   !> Header of a state: its grid layout, sites and solute coordinates, without fields
   !>
   !> @param[in] state  state
   function header_of(state) result(header)
      !> State
      type(moz_checkpoint_state_type), intent(in) :: state
      !> Header
      type(moz_checkpoint_state_type) :: header
      ! Component by component: gfortran 13 drops deferred-length strings in structure constructors
      header%kind = state%kind
      header%grid = state%grid
      header%solvent_sites = state%solvent_sites
      header%solute_sites = state%solute_sites
      if (allocated(state%solute_xyz)) header%solute_xyz = state%solute_xyz
   end function header_of

   !> Write the datasets a load reads of a Cartesian state, bypassing the checks of a save
   !>
   !> @param[in]  file   container file, replaced
   !> @param[in]  state  Cartesian state with fields
   !> @param[out] err    write failure
   subroutine write_bare_state(file, state, err)
      !> Container file
      character(len=*), intent(in) :: file
      !> State
      type(moz_checkpoint_state_type), intent(in) :: state
      !> Write failure
      type(moist_error), allocatable, intent(out) :: err
      type(binary_file) :: binf

      call delete_file(file, err)
      if (allocated(err)) return
      call binf%open(file, status="new", action="w", error=err)
      if (allocated(err)) return
      call binf%add("/checkpoint/version", 2, error=err)
      if (.not. allocated(err)) call binf%add("/checkpoint/kind", state%kind, error=err)
      if (.not. allocated(err)) call binf%add("/grid/kind", state%grid%kind, error=err)
      if (.not. allocated(err)) call binf%add("/grid/shape", state%grid%shape, error=err)
      if (.not. allocated(err)) call binf%add("/grid/spacing", state%grid%spacing, error=err)
      if (.not. allocated(err)) call binf%add("/grid/origin", state%grid%origin, error=err)
      if (.not. allocated(err)) call binf%add("/sites/solvent", state%solvent_sites, error=err)
      if (.not. allocated(err)) call binf%add("/sites/solute", state%solute_sites, error=err)
      if (.not. allocated(err)) call binf%add("/solute/coordinates", state%solute_xyz, error=err)
      if (.not. allocated(err)) call binf%add("/metadata/temperature", state%temperature, error=err)
      if (.not. allocated(err)) call binf%add("/metadata/residual", state%residual, error=err)
      if (.not. allocated(err)) call binf%add("/metadata/lambda", state%lambda, error=err)
      if (.not. allocated(err)) call binf%add("/metadata/iterations", state%iterations, error=err)
      if (.not. allocated(err)) call binf%add("/metadata/converged", merge(1, 0, state%converged), error=err)
      if (.not. allocated(err)) call binf%add("/tr", state%tr, error=err)
      if (allocated(err)) then
         call binf%close()
      else
         call binf%close(error=err)
      end if
   end subroutine write_bare_state

   !> Append a state to the shared log
   !>
   !> @param[inout] self   store
   !> @param[in]    state  state with fields
   !> @param[out]   error  store without a log
   subroutine memory_save(self, state, error)
      !> Store
      class(memory_store_type), intent(inout) :: self
      !> State
      type(moz_checkpoint_state_type), intent(in) :: state
      !> Store without a log
      type(moist_error), allocatable, intent(out) :: error
      type(moz_checkpoint_state_type), allocatable :: saved(:)
      integer :: n
      if (.not. associated(self%log)) then
         call fatal_error(error, "memory store has no log")
         return
      end if
      ! Appended without an array constructor, on which ifort 2021 aborts for these states
      n = size(self%log%saved)
      allocate (saved(n + 1))
      saved(:n) = self%log%saved
      saved(n + 1) = state
      call move_alloc(saved, self%log%saved)
   end subroutine memory_save

   !> Offer the preset state when it fits the header
   !>
   !> @param[inout] self     store
   !> @param[out]   state    preset state
   !> @param[out]   found    whether the preset fits
   !> @param[out]   source   "memory"
   !> @param[out]   warning  why the preset does not fit
   subroutine memory_load(self, state, found, source, warning)
      !> Store
      class(memory_store_type), intent(inout) :: self
      !> Preset state
      type(moz_checkpoint_state_type), intent(out) :: state
      !> Whether the preset fits
      logical, intent(out) :: found
      !> Source
      character(len=:), allocatable, intent(out) :: source
      !> Why the preset does not fit
      type(moist_error), allocatable, intent(out) :: warning

      found = .false.
      source = "memory"
      if (.not. allocated(self%preset%tr)) return
      call self%preset%check_compatible(self%header, warning)
      if (allocated(warning)) return
      state = self%preset
      found = .true.
   end subroutine memory_load

   !> Location for messages
   !>
   !> @param[in] self  store
   function memory_describe(self) result(location)
      !> Store
      class(memory_store_type), intent(in) :: self
      !> Location
      character(len=:), allocatable :: location
      location = "memory"
      if (.not. associated(self%log)) location = ""
   end function memory_describe

   !> Built-in Lebedev rule with its points in reverse order
   !>
   !> @param[in]  self     generator
   !> @param[in]  request  angular rule request
   !> @param[out] grid     reversed rule
   !> @param[out] error    invalid request
   subroutine permuted_angular_select(self, request, grid, error)
      !> Generator
      class(permuted_angular_type), intent(in) :: self
      !> Angular rule request
      type(moist_math_grid_angular_request_type), intent(in) :: request
      !> Reversed rule
      type(moist_math_grid_angular_type), intent(out) :: grid
      !> Invalid request
      type(moist_error), allocatable, intent(out) :: error

      call self%moist_math_grid_angular_generator_lebedev_type%select(request, grid, error)
      if (allocated(error)) return
      grid%points = grid%points(:, grid%npts:1:-1)
      grid%weights = grid%weights(grid%npts:1:-1)
   end subroutine permuted_angular_select

end module test_moz_checkpoint
