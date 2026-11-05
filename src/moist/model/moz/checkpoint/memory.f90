!> Checkpoint store keeping the last state in a variable of the caller
!>
!> - The store points to a state owned by the caller (the slot); models copy the
!>   store they are given, and every copy saves into and loads from the same slot
!> - Each save replaces the slot; a load offers the slot when it fits the header
!> - The caller reads the last saved solution from the slot and can fill the slot
!>   before a solve to seed it; a slot without t(r) offers nothing
!> - A field accepted without iterations is not saved, so the slot keeps the origin
!>   and solute coordinates of the solve that saved it
!> - One slot per solve family: a VV and a UV store need separate slots
!>
!>    type(moz_checkpoint_state_type), target :: slot
!>    type(moz_memory_checkpoint_type) :: store
!>
!>    call new_memory_checkpoint(store, slot)
!>    call model%set_checkpoint(store)
module moist_model_moz_checkpoint_memory
   use mctc_env, only: error_type, fatal_error
   use moist_model_moz_checkpoint, only: moz_checkpoint_type, moz_checkpoint_state_type
   implicit none(type, external)
   private

   public :: moz_memory_checkpoint_type, new_memory_checkpoint

   !> Store saving to and loading from one state of the caller
   type, extends(moz_checkpoint_type) :: moz_memory_checkpoint_type
      !> State of the caller replaced by each save; unassociated saves nothing
      type(moz_checkpoint_state_type), pointer :: slot => null()
   contains
      procedure :: save => memory_save
      procedure :: load => memory_load
      procedure :: describe => memory_describe
   end type moz_memory_checkpoint_type

contains

   !> Create a memory store on a state of the caller
   !>
   !> - The slot keeps its contents, so a state placed there before seeds the first solve
   !> - The slot must outlive every copy of the store
   !>
   !> @param[out]   self      store
   !> @param[inout] slot      state replaced by each save and offered by each load
   !> @param[in]    interval  iterations between intermediate saves; default 0
   subroutine new_memory_checkpoint(self, slot, interval)
      !> Store
      type(moz_memory_checkpoint_type), intent(out) :: self
      !> State of the caller
      type(moz_checkpoint_state_type), intent(inout), target :: slot
      !> Iterations between intermediate saves
      integer, intent(in), optional :: interval

      self%slot => slot
      if (present(interval)) self%interval = interval
   end subroutine new_memory_checkpoint

   !> Replace the slot with one state
   !>
   !> @param[inout] self   store
   !> @param[in]    state  state with fields
   !> @param[out]   error  incomplete or inconsistent state
   subroutine memory_save(self, state, error)
      !> Store
      class(moz_memory_checkpoint_type), intent(inout) :: self
      !> State with fields
      type(moz_checkpoint_state_type), intent(in) :: state
      !> Incomplete state
      type(error_type), allocatable, intent(out) :: error

      if (.not. associated(self%slot)) return
      call state%check_complete(error)
      if (allocated(error)) return
      self%slot = state
   end subroutine memory_save

   !> Offer the slot when it holds t(r) and fits the header
   !>
   !> - The field shape is taken from t(r), so a slot filled by hand needs no field_shape
   !>
   !> @param[inout] self     store with a header
   !> @param[out]   state    copy of the slot; defined when found
   !> @param[out]   found    whether the slot fits
   !> @param[out]   source   "memory"
   !> @param[out]   warning  why the slot does not fit
   subroutine memory_load(self, state, found, source, warning)
      !> Store
      class(moz_memory_checkpoint_type), intent(inout) :: self
      !> Copy of the slot
      type(moz_checkpoint_state_type), intent(out) :: state
      !> Whether the slot fits
      logical, intent(out) :: found
      !> Source
      character(len=:), allocatable, intent(out) :: source
      !> Why the slot does not fit
      type(error_type), allocatable, intent(out) :: warning
      type(error_type), allocatable :: error

      found = .false.
      source = "memory"
      if (.not. associated(self%slot)) return
      if (.not. allocated(self%slot%tr)) return
      state = self%slot
      state%field_shape = shape(state%tr)
      call state%check_compatible(self%header, error)
      if (allocated(error)) then
         call fatal_error(warning, "ignored memory ("//error%message//")")
         return
      end if
      found = .true.
   end subroutine memory_load

   !> Location for messages
   !>
   !> @param[in] self  store
   function memory_describe(self) result(location)
      !> Store
      class(moz_memory_checkpoint_type), intent(in) :: self
      !> "memory"; empty when nothing is saved
      character(len=:), allocatable :: location
      location = ""
      if (associated(self%slot)) location = "memory"
   end function memory_describe

end module moist_model_moz_checkpoint_memory
