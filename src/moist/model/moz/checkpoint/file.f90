!> Checkpoint store backed by binary container files (.mbin)
!>
!> - Saves write a temporary file `<file>.<n>.tmp` of their own and rename it over
!>   `file`, so an interrupted or failed save keeps the previous checkpoint and
!>   concurrent saves never write into each other's file; an empty name saves nothing
!> - Loads try `guesses` in order and take the first existing, readable file
!>   that fits the header; missing files are skipped silently
!> - Loads check the metadata and the extents of /tr against the header before
!>   reading /tr
!> - Layout (format version 2), lengths in bohr and temperatures in K:
!>   - /checkpoint/version, /checkpoint/kind
!>   - /grid/kind, /grid/shape, /grid/spacing; Cartesian /grid/origin; radial /grid/r,
!>     /grid/weights; point grids /grid/xyz, /grid/weights, /grid/owner and the
!>     Fourier basis /grid/kspace/modes, /grid/kspace/spacing, /grid/kspace/ref
!>   - atom-centered grids: /grid/points, /grid/shells/offset, /grid/shells/radius,
!>     /grid/shells/points, /grid/angular
!>   - /sites/solvent, optional /sites/solute and /solute/coordinates
!>   - 3D UV: /longrange/status ("known", "zero" or "unknown"); a known tail also
!>     /longrange/alpha, /charges/solute, /charges/solvent
!>   - /metadata/closure and /metadata/theory (both optional), /metadata/temperature,
!>     /metadata/residual, /metadata/lambda, /metadata/iterations, /metadata/converged
!>   - /tr: physical indirect correlation, (npoint, ncol)
module moist_model_moz_checkpoint_file
   use mctc_env, only: wp, error_type, fatal_error
   use moist_utils_binio, only: binary_file
   use moist_utils_file, only: create_new_file, replace_file, delete_file
   use moist_model_moz_checkpoint, only: moz_checkpoint_type, moz_checkpoint_state_type
   implicit none(type, external)
   private

   public :: moz_file_checkpoint_type, new_file_checkpoint

   !> Format version written into /checkpoint/version
   integer, parameter :: format_version = 2

   !> Store saving to one file and loading from candidate files
   type, extends(moz_checkpoint_type) :: moz_file_checkpoint_type
      !> File replaced by each save; empty saves nothing
      character(len=:), allocatable :: file
      !> Files tried by load, in order; blank entries are skipped
      character(len=:), allocatable :: guesses(:)
   contains
      procedure :: save => file_save
      procedure :: load => file_load
      procedure :: describe => file_describe
   end type moz_file_checkpoint_type

contains

   !> Create a file store
   !>
   !> @param[out] self      store
   !> @param[in]  file      file replaced by each save; empty saves nothing
   !> @param[in]  guesses   files tried by load, in order; default `file`
   !> @param[in]  interval  iterations between intermediate saves; default 0
   subroutine new_file_checkpoint(self, file, guesses, interval)
      !> Store
      type(moz_file_checkpoint_type), intent(out) :: self
      !> Saved file
      character(len=*), intent(in) :: file
      !> Loaded files
      character(len=*), intent(in), optional :: guesses(:)
      !> Iterations between intermediate saves
      integer, intent(in), optional :: interval

      self%file = trim(file)
      if (present(guesses)) then
         allocate (character(len=len(guesses)) :: self%guesses(size(guesses)))
         self%guesses(:) = guesses
      else
         self%guesses = [character(len=len(self%file)) :: self%file]
      end if
      if (present(interval)) self%interval = interval
   end subroutine new_file_checkpoint

   !> Replace the file with one state, keeping the previous file on failure
   !>
   !> @param[inout] self   store
   !> @param[in]    state  state with fields
   !> @param[out]   error  incomplete or inconsistent state, no temporary file, or open, write,
   !>                       close or rename failure; a failed cleanup of the temporary file is
   !>                       added to the message
   subroutine file_save(self, state, error)
      !> Store
      class(moz_file_checkpoint_type), intent(inout) :: self
      !> State with fields
      type(moz_checkpoint_state_type), intent(in) :: state
      !> Failure
      type(error_type), allocatable, intent(out) :: error
      type(binary_file) :: binf
      type(error_type), allocatable :: cleanup
      character(len=:), allocatable :: partial

      if (.not. allocated(self%file)) return
      if (len(self%file) == 0) return
      call state%check_complete(error)
      if (allocated(error)) return
      call create_new_file(self%file, partial, error)
      if (allocated(error)) return
      call binf%open(partial, status="old", action="w", error=error, initialize=.true.)
      if (allocated(error)) then
         call delete_file(partial, cleanup)
         return
      end if
      call write_state(binf, state, error)
      if (allocated(error)) then
         call binf%close()
      else
         call binf%close(error=error)
      end if
      if (.not. allocated(error)) call replace_file(partial, self%file, error)
      if (allocated(error)) then
         call delete_file(partial, cleanup)
         if (allocated(cleanup)) error%message = error%message//"; "//cleanup%message
      end if
   end subroutine file_save

   !> Load the first existing, readable candidate fitting the header
   !>
   !> @param[inout] self     store with a header
   !> @param[out]   state    loaded state; defined when found
   !> @param[out]   found    whether a fitting state was loaded
   !> @param[out]   source   file of the loaded state
   !> @param[out]   warning  existing candidates that were unreadable or did not fit
   subroutine file_load(self, state, found, source, warning)
      !> Store
      class(moz_file_checkpoint_type), intent(inout) :: self
      !> Loaded state
      type(moz_checkpoint_state_type), intent(out) :: state
      !> Whether a state was loaded
      logical, intent(out) :: found
      !> File of the loaded state
      character(len=:), allocatable, intent(out) :: source
      !> Skipped candidates
      type(error_type), allocatable, intent(out) :: warning
      type(error_type), allocatable :: error
      character(len=:), allocatable :: skipped, file
      logical :: exists
      integer :: i

      found = .false.
      source = ""
      skipped = ""
      if (.not. allocated(self%guesses)) return
      do i = 1, size(self%guesses)
         file = trim(self%guesses(i))
         if (len(file) == 0) cycle
         inquire (file=file, exist=exists)
         if (.not. exists) cycle
         call read_state(file, self%header, state, error)
         if (allocated(error)) then
            if (len(skipped) > 0) skipped = skipped//"; "
            skipped = skipped//file//" ("//error%message//")"
            deallocate (error)
            cycle
         end if
         found = .true.
         source = file
         exit
      end do
      if (len(skipped) > 0) call fatal_error(warning, "ignored "//skipped)
   end subroutine file_load

   !> File replaced by each save
   !>
   !> @param[in] self  store
   function file_describe(self) result(location)
      !> Store
      class(moz_file_checkpoint_type), intent(in) :: self
      !> Saved file; empty when nothing is saved
      character(len=:), allocatable :: location
      location = ""
      if (allocated(self%file)) location = self%file
   end function file_describe

   !> Write the datasets of one state into an open container
   !>
   !> @param[inout] binf   container open for writing
   !> @param[in]    state  state with fields
   !> @param[out]   error  first write failure
   subroutine write_state(binf, state, error)
      !> Container
      type(binary_file), intent(inout) :: binf
      !> State
      type(moz_checkpoint_state_type), intent(in) :: state
      !> First failure
      type(error_type), allocatable, intent(out) :: error

      call binf%add("/checkpoint/version", format_version, error=error)
      if (.not. allocated(error)) call binf%add("/checkpoint/kind", state%kind, error=error)
      if (allocated(state%grid%kind) .and. .not. allocated(error)) then
         call binf%add("/grid/kind", state%grid%kind, error=error)
      end if
      if (.not. allocated(error)) call binf%add("/grid/shape", state%grid%shape, error=error)
      if (.not. allocated(error)) call binf%add("/grid/spacing", state%grid%spacing, error=error)
      if (.not. allocated(error)) call binf%adda("/grid/spacing", "units", "bohr", error=error)
      if (allocated(state%grid%origin) .and. .not. allocated(error)) then
         call binf%add("/grid/origin", state%grid%origin, error=error)
         if (.not. allocated(error)) call binf%adda("/grid/origin", "units", "bohr", error=error)
      end if
      if (.not. allocated(error)) call write_points(binf, state, error)
      if (allocated(state%grid%angular_grid) .and. .not. allocated(error)) then
         call binf%add("/grid/angular", state%grid%angular_grid, error=error)
      end if
      if (allocated(state%grid%point_index) .and. .not. allocated(error)) then
         call binf%add("/grid/points", state%grid%point_index, error=error)
      end if
      if (allocated(state%grid%shell_offset) .and. allocated(state%grid%shell_radius) .and. allocated(state%grid%shell_points) &
         & .and. .not. allocated(error)) then
         call binf%add("/grid/shells/offset", state%grid%shell_offset, error=error)
         if (.not. allocated(error)) call binf%add("/grid/shells/radius", state%grid%shell_radius, error=error)
         if (.not. allocated(error)) call binf%adda("/grid/shells/radius", "units", "bohr", error=error)
         if (.not. allocated(error)) call binf%add("/grid/shells/points", state%grid%shell_points, error=error)
      end if
      if (.not. allocated(error)) call binf%add("/sites/solvent", state%solvent_sites, error=error)
      if (allocated(state%solute_sites) .and. .not. allocated(error)) then
         call binf%add("/sites/solute", state%solute_sites, error=error)
      end if
      if (allocated(state%solute_xyz) .and. .not. allocated(error)) then
         call binf%add("/solute/coordinates", state%solute_xyz, error=error)
         if (.not. allocated(error)) call binf%adda("/solute/coordinates", "units", "bohr", error=error)
      end if
      if (allocated(state%tail) .and. .not. allocated(error)) then
         call binf%add("/longrange/status", state%tail, error=error)
         if (allocated(state%solute_charges) .and. allocated(state%solvent_charges) .and. .not. allocated(error)) then
            call binf%add("/longrange/alpha", state%alpha, error=error)
            if (.not. allocated(error)) call binf%adda("/longrange/alpha", "units", "1/bohr", error=error)
            if (.not. allocated(error)) call binf%add("/charges/solute", state%solute_charges, error=error)
            if (.not. allocated(error)) call binf%add("/charges/solvent", state%solvent_charges, error=error)
         end if
      end if
      if (allocated(state%closure) .and. .not. allocated(error)) then
         call binf%add("/metadata/closure", state%closure, error=error)
      end if
      if (allocated(state%theory) .and. .not. allocated(error)) then
         call binf%add("/metadata/theory", state%theory, error=error)
      end if
      if (.not. allocated(error)) call binf%add("/metadata/temperature", state%temperature, error=error)
      if (.not. allocated(error)) call binf%adda("/metadata/temperature", "units", "K", error=error)
      if (.not. allocated(error)) call binf%add("/metadata/residual", state%residual, error=error)
      if (.not. allocated(error)) call binf%add("/metadata/lambda", state%lambda, error=error)
      if (.not. allocated(error)) call binf%add("/metadata/iterations", state%iterations, error=error)
      if (.not. allocated(error)) call binf%add("/metadata/converged", merge(1, 0, state%converged), error=error)
      if (.not. allocated(error)) call binf%add("/tr", state%tr, error=error)
      if (.not. allocated(error)) call binf%adda("/tr", "description", "indirect_correlation_function", error=error)
   end subroutine write_state

   !> Write the nodes, points, weights, owners and Fourier basis of a layout
   !>
   !> @param[inout] binf   container open for writing
   !> @param[in]    state  state with its layout
   !> @param[out]   error  first write failure
   subroutine write_points(binf, state, error)
      !> Container
      type(binary_file), intent(inout) :: binf
      !> State
      type(moz_checkpoint_state_type), intent(in) :: state
      !> First failure
      type(error_type), allocatable, intent(out) :: error

      associate (grid => state%grid)
         if (allocated(grid%r)) then
            call binf%add("/grid/r", grid%r, error=error)
            if (.not. allocated(error)) call binf%adda("/grid/r", "units", "bohr", error=error)
         end if
         if (allocated(grid%xyz) .and. .not. allocated(error)) then
            call binf%add("/grid/xyz", grid%xyz, error=error)
            if (.not. allocated(error)) call binf%adda("/grid/xyz", "units", "bohr", error=error)
         end if
         if (allocated(grid%w) .and. .not. allocated(error)) call binf%add("/grid/weights", grid%w, error=error)
         if (allocated(grid%owner) .and. .not. allocated(error)) call binf%add("/grid/owner", grid%owner, error=error)
         if (allocated(grid%kmodes) .and. allocated(grid%kspacing) .and. allocated(grid%kref) &
            & .and. .not. allocated(error)) then
            call binf%add("/grid/kspace/modes", grid%kmodes, error=error)
            if (.not. allocated(error)) call binf%add("/grid/kspace/spacing", grid%kspacing, error=error)
            if (.not. allocated(error)) call binf%adda("/grid/kspace/spacing", "units", "1/bohr", error=error)
            if (.not. allocated(error)) call binf%add("/grid/kspace/ref", grid%kref, error=error)
            if (.not. allocated(error)) call binf%adda("/grid/kspace/ref", "units", "bohr", error=error)
         end if
      end associate
   end subroutine write_points

   !> Read one state fitting a header from a file
   !>
   !> - Metadata and the extents of /tr are checked before /tr is read
   !>
   !> @param[in]  file    container file
   !> @param[in]  header  discretization of the current solve
   !> @param[out] state   state with t(r)
   !> @param[out] error   open, read, version or compatibility failure
   subroutine read_state(file, header, state, error)
      !> Container file
      character(len=*), intent(in) :: file
      !> Discretization of the current solve
      type(moz_checkpoint_state_type), intent(in) :: header
      !> State
      type(moz_checkpoint_state_type), intent(out) :: state
      !> Failure
      type(error_type), allocatable, intent(out) :: error
      type(binary_file) :: binf

      call binf%open(file, status="old", action="r", error=error)
      if (allocated(error)) return
      call read_metadata(binf, state, error)
      if (.not. allocated(error)) call state%check_compatible(header, error)
      if (.not. allocated(error)) call binf%get("/tr", state%tr, error=error)
      if (allocated(error)) then
         call binf%close()
         return
      end if
      call binf%close(error=error)
   end subroutine read_state

   !> Read the metadata of one state and the extents of its fields
   !>
   !> @param[inout] binf   container open for reading
   !> @param[inout] state  state receiving the metadata
   !> @param[out]   error  first read or version failure
   subroutine read_metadata(binf, state, error)
      !> Container
      type(binary_file), intent(inout) :: binf
      !> State
      type(moz_checkpoint_state_type), intent(inout) :: state
      !> First failure
      type(error_type), allocatable, intent(out) :: error
      integer :: version, converged

      call binf%get("/checkpoint/version", version, error=error)
      if (.not. allocated(error)) then
         if (version /= format_version) call fatal_error(error, "unsupported checkpoint format version")
      end if
      if (.not. allocated(error)) call binf%get("/checkpoint/kind", state%kind, error=error)
      if (.not. allocated(error)) call binf%get("/grid/shape", state%grid%shape, error=error)
      if (.not. allocated(error)) call binf%get("/grid/spacing", state%grid%spacing, error=error)
      if (.not. allocated(error)) call binf%get("/sites/solvent", state%solvent_sites, error=error)
      if (.not. allocated(error)) call binf%get_shape("/tr", state%field_shape, error=error)
      if (.not. allocated(error)) call binf%get("/metadata/temperature", state%temperature, error=error)
      if (.not. allocated(error)) call binf%get("/metadata/residual", state%residual, error=error)
      if (.not. allocated(error)) call binf%get("/metadata/lambda", state%lambda, error=error)
      if (.not. allocated(error)) call binf%get("/metadata/iterations", state%iterations, error=error)
      if (.not. allocated(error)) then
         call binf%get("/metadata/converged", converged, error=error)
         state%converged = converged /= 0
      end if
      if (.not. allocated(error)) call read_optional(binf, state, error)
   end subroutine read_metadata

   !> Read the optional metadata of one state
   !>
   !> @param[inout] binf   container open for reading
   !> @param[inout] state  state receiving the optional fields
   !> @param[out]   error  first read failure
   subroutine read_optional(binf, state, error)
      !> Container
      type(binary_file), intent(inout) :: binf
      !> State
      type(moz_checkpoint_state_type), intent(inout) :: state
      !> First failure
      type(error_type), allocatable, intent(out) :: error

      if (has(binf, "/grid/kind", error)) call binf%get("/grid/kind", state%grid%kind, error=error)
      if (has(binf, "/grid/origin", error)) call binf%get("/grid/origin", state%grid%origin, error=error)
      if (has(binf, "/grid/r", error)) call binf%get("/grid/r", state%grid%r, error=error)
      if (has(binf, "/grid/xyz", error)) call binf%get("/grid/xyz", state%grid%xyz, error=error)
      if (has(binf, "/grid/weights", error)) call binf%get("/grid/weights", state%grid%w, error=error)
      if (has(binf, "/grid/owner", error)) call binf%get("/grid/owner", state%grid%owner, error=error)
      if (has(binf, "/grid/kspace/modes", error)) then
         call binf%get("/grid/kspace/modes", state%grid%kmodes, error=error)
         if (.not. allocated(error)) call binf%get("/grid/kspace/spacing", state%grid%kspacing, error=error)
         if (.not. allocated(error)) call binf%get("/grid/kspace/ref", state%grid%kref, error=error)
      end if
      if (has(binf, "/grid/angular", error)) call binf%get("/grid/angular", state%grid%angular_grid, error=error)
      if (has(binf, "/grid/points", error)) call binf%get("/grid/points", state%grid%point_index, error=error)
      if (has(binf, "/grid/shells/offset", error)) then
         call binf%get("/grid/shells/offset", state%grid%shell_offset, error=error)
         if (.not. allocated(error)) call binf%get("/grid/shells/radius", state%grid%shell_radius, error=error)
         if (.not. allocated(error)) call binf%get("/grid/shells/points", state%grid%shell_points, error=error)
      end if
      if (has(binf, "/sites/solute", error)) call binf%get("/sites/solute", state%solute_sites, error=error)
      if (has(binf, "/solute/coordinates", error)) call binf%get("/solute/coordinates", state%solute_xyz, error=error)
      if (has(binf, "/metadata/closure", error)) call binf%get("/metadata/closure", state%closure, error=error)
      if (has(binf, "/metadata/theory", error)) call binf%get("/metadata/theory", state%theory, error=error)
      if (has(binf, "/longrange/status", error)) call binf%get("/longrange/status", state%tail, error=error)
      if (has(binf, "/longrange/alpha", error)) then
         call binf%get("/longrange/alpha", state%alpha, error=error)
         if (.not. allocated(error)) call binf%get("/charges/solute", state%solute_charges, error=error)
         if (.not. allocated(error)) call binf%get("/charges/solvent", state%solvent_charges, error=error)
      end if
   end subroutine read_optional

   !> Whether a container holds a path; false once an earlier read failed
   !>
   !> @param[inout] binf   container open for reading
   !> @param[in]    path   dataset path
   !> @param[inout] error  earlier or lookup failure
   function has(binf, path, error) result(exists)
      !> Container
      type(binary_file), intent(inout) :: binf
      !> Dataset path
      character(len=*), intent(in) :: path
      !> Earlier or lookup failure
      type(error_type), allocatable, intent(inout) :: error
      !> Whether the path exists and no failure occurred
      logical :: exists

      exists = .false.
      if (allocated(error)) return
      call binf%exist(path, exists, error=error)
      exists = exists .and. .not. allocated(error)
   end function has

end module moist_model_moz_checkpoint_file
