!> Tests of the MOZ checkpoints: states, the file store and warm starts: moz_checkpoint
!>
!> - File store round trip, candidate order, refusals and failed saves on small synthetic states
!> - cSPC/E water VV and water-in-water 1D and 3D UV solves with a memory store;
!>   a solve that iterates saves its final field, so a reload without a new save
!>   was accepted without iterations
!> - Each test runs as a selected test
module test_moz_checkpoint
   use mctc_env, only: wp, moist_error => error_type, fatal_error
   use mctc_io, only: structure_type
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_context, only: moist_context_type, new_context
   use moist_model_type, only: solvation_model_type
   use moist_channels_coupling, only: coupling_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_model_moz_checkpoint, only: moz_checkpoint_type, moz_checkpoint_state_type, &
      & checkpoint_kind_vv, checkpoint_kind_uv1d, checkpoint_kind_uv3d
   use moist_model_moz_checkpoint_file, only: moz_file_checkpoint_type, new_file_checkpoint
   use moist_model_moz_solvent_vv, only: solvent_vv_type, solvent_vv_parameters_type, solvent_vv_thermodynamics_type
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_potential, only: lj_cspce, coulomb_spce
   use test_helpers, only: check_moist_error
   use test_moz_fixtures, only: new_cspce_water_solvent
   implicit none(type, external)
   private
   public :: collect_moz_checkpoint

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
      & file_trip = "moist-test-checkpoint-trip.mbin", file_save = "moist-test-checkpoint-save.mbin"

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
         & new_unittest("vv_reload", check_vv_reload), &
         & new_unittest("vv_intermediate", check_vv_intermediate), &
         & new_unittest("vv_stages", check_vv_stages), &
         & new_unittest("vv_unfitting", check_vv_unfitting), &
         & new_unittest("uv1d_reload", check_uv1d_reload), &
         & new_unittest("uv3d_reload", check_uv3d_reload)]
   end subroutine collect_moz_checkpoint

   !> Save a state to a file and read it back unchanged
   !>
   !> - Every dataset survives, including the optional provenance
   !> - An empty file name saves nothing; describe names the saved file
   subroutine check_file_round_trip(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
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
      call delete_file(file_trip)
      call check(error, found .and. .not. allocated(warning) .and. source == file_trip, more="saved state not found")
      if (allocated(error)) return
      call check(error, loaded%kind == state%kind .and. all(loaded%grid_shape == state%grid_shape) &
         & .and. loaded%spacing == state%spacing .and. all(loaded%origin == state%origin) &
         & .and. all(loaded%solvent_sites == state%solvent_sites) .and. all(loaded%solute_sites == state%solute_sites) &
         & .and. all(loaded%solute_xyz == state%solute_xyz), more="discretization or sites changed")
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
   !> - Unreadable files, states of another grid or sites, and fields of another
   !>   shape than the header expects are skipped with a warning
   !> - The first fitting file is loaded, also as an initial guess
   subroutine check_file_candidates(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
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
      other%grid_shape = [5]
      other%tr = other%tr(:5, :)
      malformed = state
      malformed%tr = malformed%tr(:7, :)
      call new_file_checkpoint(writer, file_b)
      call writer%save(other, err)
      if (.not. allocated(err)) then
         call new_file_checkpoint(writer, file_d)
         call writer%save(malformed, err)
      end if
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
      call check(error, all(loaded%tr == state%tr), more="loaded the wrong fields")
      if (allocated(error)) return

      store%header = header_of(state)
      call store%load_guess(ctx, 1.0_wp, spread([0.0_wp, 0.0_wp], 1, 8), 4, guess)
      call delete_file(file_a)
      call delete_file(file_b)
      call delete_file(file_d)
      call check(error, allocated(guess), more="the malformed candidate hid the fitting one from the guess")
      if (allocated(error)) return
      call check(error, all(guess == state%tr), more="the guess holds the wrong fields")
      if (allocated(error)) return

      store%header%solvent_sites = [8, 1, 1]
      call store%load(loaded, found, source, warning)
      call check(error, .not. found .and. allocated(warning), more="a state of other sites was loaded")
      if (allocated(error)) return

      call delete_file(file_c)
      store%header = header_of(state)
      call store%load(loaded, found, source, warning)
      call check(error, .not. found .and. .not. allocated(warning), more="missing files were not skipped silently")
   end subroutine check_file_candidates

   !> Saves replace the file whole and keep it on failure
   !>
   !> - A stale temporary file of an interrupted save does not block a save
   !> - A later save replaces the existing file
   !> - An incomplete state is refused and the previous file stays loadable
   !> - No temporary file is left behind
   subroutine check_file_failed_save(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err, warning
      type(moz_file_checkpoint_type) :: store
      type(moz_checkpoint_state_type) :: state, incomplete, loaded
      character(len=:), allocatable :: source
      logical :: found, partial
      integer :: unit

      call synthetic_state(state)
      call new_file_checkpoint(store, file_save)
      store%header = header_of(state)
      open (newunit=unit, file=file_save//".tmp", status="replace", action="write")
      write (unit, "(a)") "interrupted save"
      close (unit)
      call store%save(state, err)
      if (.not. allocated(err)) then
         state%residual = 2.5e-10_wp
         call store%save(state, err)
      end if
      call check_moist_error(error, err)
      if (allocated(error)) then
         call delete_file(file_save)
         call delete_file(file_save//".tmp")
         return
      end if

      incomplete = state
      deallocate (incomplete%tr)
      incomplete%residual = 1.0_wp
      call store%save(incomplete, err)
      call store%load(loaded, found, source, warning)
      inquire (file=file_save//".tmp", exist=partial)
      call delete_file(file_save)
      call delete_file(file_save//".tmp")
      call check(error, allocated(err), more="an incomplete state was saved")
      if (allocated(error)) return
      call check(error, found .and. .not. allocated(warning), more="the failed save lost the previous file")
      if (allocated(error)) return
      call check(error, loaded%residual == state%residual, more="the later save did not replace the file")
      if (allocated(error)) return
      call check(error, .not. partial, more="a temporary file was left behind")
   end subroutine check_file_failed_save

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
            & .and. all(saved%grid_shape == [params%npts]) .and. all(saved%solvent_sites == [8, 1, 1]) &
            & .and. saved%closure == "KH", more="unexpected saved VV state")
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

      store%preset = log%saved(n - 2)
      deallocate (log%saved)
      allocate (log%saved(0))
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

   !> Coupling stages: saves of the last stage only, and a guess skips the ramp
   !>
   !> - Two stages save only fields at full coupling
   !> - A converged guess is accepted without iterations despite the stages
   subroutine check_vv_stages(error)
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
      params%nlambda = 2
      call new_memory_store(store, log)
      store%interval = 1
      call solve_vv(ctx, params, store, cold, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      n = size(log%saved)
      call check(error, n >= 1, more="the staged solve saved nothing")
      if (allocated(error)) return
      call check(error, all(log%saved(:)%lambda == 1.0_wp) .and. log%saved(n)%converged, &
         & more="a state below full coupling was saved")
      if (allocated(error)) return

      store%preset = log%saved(n)
      call solve_vv(ctx, params, store, warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == n .and. warm%is_solved(params), &
         & more="a converged guess ran the coupling stages")
   end subroutine check_vv_stages

   !> States of another grid or kind are ignored; the solve starts cold
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
      store%preset%grid_shape = [params%npts/2]
      call solve_vv(ctx, params, store, warm, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, size(log%saved) == 2 .and. log%saved(2)%iterations > 0 .and. warm%is_solved(params), &
         & more="a state of another grid was used")
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
            & .and. all(saved%solute_sites == [8, 1, 1]) .and. all(saved%grid_shape == [params%solvent%npts]), &
            & more="unexpected saved 1D UV state")
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
            & .and. all(saved%grid_shape == [32, 32, 32]) .and. allocated(saved%origin), &
            & more="unexpected saved 3D UV state")
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

   !> Small grids and tolerance of the VV solves
   !>
   !> @param[out] params  settings
   subroutine vv_settings(params)
      !> Settings
      type(solvent_vv_parameters_type), intent(out) :: params
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
      type(moist_math_grid_3d_cartesian_type), intent(in) :: grid
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

   !> Synthetic 3D UV state with every optional dataset
   !>
   !> @param[out] state  state with fields (8, 2)
   subroutine synthetic_state(state)
      !> State
      type(moz_checkpoint_state_type), intent(out) :: state
      integer :: i, j

      state%kind = checkpoint_kind_uv3d
      state%grid_shape = [2, 2, 2]
      state%spacing = 0.5_wp
      state%origin = [-1.0_wp, -2.0_wp, -3.0_wp]
      state%solvent_sites = [8, 1]
      state%solute_sites = [17]
      state%solute_xyz = reshape([0.1_wp, 0.2_wp, 0.3_wp], [3, 1])
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

   !> Header of a state: its discretization and sites, without fields
   !>
   !> @param[in] state  state
   function header_of(state) result(header)
      !> State
      type(moz_checkpoint_state_type), intent(in) :: state
      !> Header
      type(moz_checkpoint_state_type) :: header
      header = moz_checkpoint_state_type(kind=state%kind, grid_shape=state%grid_shape, spacing=state%spacing, &
         & solvent_sites=state%solvent_sites, solute_sites=state%solute_sites)
   end function header_of

   !> Delete a file when it exists
   !>
   !> @param[in] file  file name
   subroutine delete_file(file)
      !> File name
      character(len=*), intent(in) :: file
      integer :: unit, stat
      logical :: exists
      inquire (file=file, exist=exists)
      if (.not. exists) return
      open (newunit=unit, file=file, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine delete_file

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
      if (.not. associated(self%log)) then
         call fatal_error(error, "memory store has no log")
         return
      end if
      self%log%saved = [self%log%saved, state]
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

end module test_moz_checkpoint
