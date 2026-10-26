!> Checkpoint store backed by binary container files (.mbin)
!>
!> - Saves write `<file>.tmp` and rename it over `file`, so an interrupted or
!>   failed save keeps the previous checkpoint; an empty name saves nothing
!> - Loads try `guesses` in order and take the first existing, readable file
!>   that fits the header; missing files are skipped silently
!> - Loads check the metadata and the extents of /tr against the header before
!>   reading /tr
!> - Layout (format version 1), lengths in bohr and temperatures in K:
!>   - /checkpoint/version, /checkpoint/kind
!>   - /grid/shape, /grid/spacing, optional /grid/origin
!>   - /sites/solvent, optional /sites/solute and /solute/coordinates
!>   - /metadata/closure, /metadata/theory (optional), /metadata/temperature,
!>     /metadata/residual, /metadata/lambda, /metadata/iterations, /metadata/converged
!>   - /tr: physical indirect correlation, (npoint, ncol)
module moist_model_moz_checkpoint_file
   use mctc_env, only: wp, error_type, fatal_error
   use moist_utils_binio, only: binary_file
   use moist_utils_file, only: replace_file, delete_file
   use moist_model_moz_checkpoint, only: moz_checkpoint_type, moz_checkpoint_state_type
   implicit none(type, external)
   private

   public :: moz_file_checkpoint_type, new_file_checkpoint

   !> Format version written into /checkpoint/version
   integer, parameter :: format_version = 1

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
   !> @param[out]   error  incomplete state, or open, write, close or rename failure
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
      if (.not. allocated(state%kind) .or. .not. allocated(state%grid_shape) .or. &
         & .not. allocated(state%solvent_sites) .or. .not. allocated(state%tr)) then
         call fatal_error(error, "Checkpoint state needs a kind, grid shape, solvent sites and t(r)")
         return
      end if
      partial = self%file//".tmp"
      call binf%open(partial, status="replace", action="w", error=error)
      if (allocated(error)) return
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
      if (.not. allocated(error)) call binf%add("/grid/shape", state%grid_shape, error=error)
      if (.not. allocated(error)) call binf%add("/grid/spacing", state%spacing, error=error)
      if (.not. allocated(error)) call binf%adda("/grid/spacing", "units", "bohr", error=error)
      if (allocated(state%origin) .and. .not. allocated(error)) then
         call binf%add("/grid/origin", state%origin, error=error)
         if (.not. allocated(error)) call binf%adda("/grid/origin", "units", "bohr", error=error)
      end if
      if (.not. allocated(error)) call binf%add("/sites/solvent", state%solvent_sites, error=error)
      if (allocated(state%solute_sites) .and. .not. allocated(error)) then
         call binf%add("/sites/solute", state%solute_sites, error=error)
      end if
      if (allocated(state%solute_xyz) .and. .not. allocated(error)) then
         call binf%add("/solute/coordinates", state%solute_xyz, error=error)
         if (.not. allocated(error)) call binf%adda("/solute/coordinates", "units", "bohr", error=error)
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
      if (.not. allocated(error)) call binf%get("/grid/shape", state%grid_shape, error=error)
      if (.not. allocated(error)) call binf%get("/grid/spacing", state%spacing, error=error)
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

      if (has(binf, "/grid/origin", error)) call binf%get("/grid/origin", state%origin, error=error)
      if (has(binf, "/sites/solute", error)) call binf%get("/sites/solute", state%solute_sites, error=error)
      if (has(binf, "/solute/coordinates", error)) call binf%get("/solute/coordinates", state%solute_xyz, error=error)
      if (has(binf, "/metadata/closure", error)) call binf%get("/metadata/closure", state%closure, error=error)
      if (has(binf, "/metadata/theory", error)) call binf%get("/metadata/theory", state%theory, error=error)
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
