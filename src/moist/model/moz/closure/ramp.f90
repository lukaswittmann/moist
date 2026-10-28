!> Cold-start closure ramp from KH up to the requested closure
!>
!> - Ladder KH, PSE2, PSE3, PSE4, HNC; PSE-n reduces to KH for n = 1 and tends to HNC
!> - A cold solve starts on KH and steps up one closure whenever the RMS
!>   residual falls below the switch threshold
!> - Closures off the ladder, or a disabled ramp, run alone
!> - Every rung is built once at setup, so switching cannot fail
module moist_model_moz_closure_ramp
   use mctc_env, only: wp, error_type
   use mctc_io_utils, only: to_lower
   use moist_context, only: moist_context_type
   use moist_model_moz_closure, only: moist_model_moz_closure_type, new_closure
   implicit none(type, external)
   private

   public :: closure_ladder, moist_model_moz_closure_ramp_type, new_closure_ramp

   !> Closures from the softest to the hardest
   character(len=4), parameter :: closure_ladder(4) = [character(len=4) :: "KH", "PSE2", "PSE3", "HNC"]

   !> Closure of one rung
   type :: closure_rung_type
      !> Closure
      class(moist_model_moz_closure_type), allocatable :: closure
   end type closure_rung_type

   !> Rungs of the ramp and the current one
   type :: moist_model_moz_closure_ramp_type
      !> Closures of the rungs; the last one is requested
      type(closure_rung_type), allocatable, private :: rungs(:)
      !> Closure names of the rungs
      character(len=8), allocatable :: names(:)
      !> Current rung
      integer :: rung = 0
      !> RMS residual below which the next rung starts
      real(wp) :: rmse = 0.0_wp
   contains
      !> Put the closure on the first rung of a cold solve
      procedure :: start => ramp_start
      !> Put the closure back on the requested rung
      procedure :: finish => ramp_finish
      !> Step up one rung once the residual is below the threshold
      procedure :: step => ramp_step
      !> Whether the current rung is the requested closure
      procedure :: at_target => ramp_at_target
      !> Name of the current closure
      procedure :: name => ramp_name
      !> Table note on the last step
      procedure :: note => ramp_note
   end type moist_model_moz_closure_ramp_type

contains

   !> Build the rungs up to the requested closure, starting on the requested one
   !>
   !> - Cites every rung on the context; a cold solve uses all of them
   !>
   !> @param[out] self    ramp on its requested rung
   !> @param[in]  name    requested closure: HNC, KH or PSE2-4
   !> @param[in]  enabled whether cold solves ramp up from KH
   !> @param[in]  rmse    RMS residual below which the next rung starts
   !> @param[out] error   unknown closure identifier
   !> @param[in,out] ctx  run context for the OpenMP team size and citations
   subroutine new_closure_ramp(self, name, enabled, rmse, error, ctx)
      !> Ramp
      type(moist_model_moz_closure_ramp_type), intent(out) :: self
      !> Requested closure
      character(len=*), intent(in) :: name
      !> Whether cold solves ramp up
      logical, intent(in) :: enabled
      !> Switch threshold
      real(wp), intent(in) :: rmse
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Run context
      type(moist_context_type), intent(inout), optional :: ctx
      integer :: i, top

      top = 0
      if (enabled) then
         do i = 1, size(closure_ladder)
            if (to_lower(trim(closure_ladder(i))) == to_lower(trim(adjustl(name)))) top = i
         end do
      end if
      if (top > 0) then
         self%names = closure_ladder(:top)
      else
         self%names = [character(len=8) :: adjustl(name)]
      end if
      allocate (self%rungs(size(self%names)))
      do i = 1, size(self%names)
         call new_closure(self%rungs(i)%closure, self%names(i), error, ctx)
         if (allocated(error)) return
      end do
      self%rung = size(self%rungs)
      self%rmse = rmse
   end subroutine new_closure_ramp

   !> Put the closure on the first rung of a cold solve
   !>
   !> @param[inout] self    ramp
   !> @param[inout] closure closure of the solve
   subroutine ramp_start(self, closure)
      !> Ramp
      class(moist_model_moz_closure_ramp_type), intent(inout) :: self
      !> Closure of the solve
      class(moist_model_moz_closure_type), allocatable, intent(inout) :: closure

      call select_rung(self, closure, 1)
   end subroutine ramp_start

   !> Put the closure back on the requested rung
   !>
   !> @param[inout] self    ramp
   !> @param[inout] closure closure of the solve
   subroutine ramp_finish(self, closure)
      !> Ramp
      class(moist_model_moz_closure_ramp_type), intent(inout) :: self
      !> Closure of the solve
      class(moist_model_moz_closure_type), allocatable, intent(inout) :: closure

      call select_rung(self, closure, size(self%rungs))
   end subroutine ramp_finish

   !> Step up one rung once the residual is below the threshold
   !>
   !> - Also steps below the solve tolerance, so a lower rung never converges
   !>
   !> @param[inout] self      ramp
   !> @param[inout] closure   closure of the solve
   !> @param[in]    rmse      RMS residual of the last iterate
   !> @param[in]    tolerance RMS residual of a converged solve
   !> @param[out]   switched  whether the closure changed
   subroutine ramp_step(self, closure, rmse, tolerance, switched)
      !> Ramp
      class(moist_model_moz_closure_ramp_type), intent(inout) :: self
      !> Closure of the solve
      class(moist_model_moz_closure_type), allocatable, intent(inout) :: closure
      !> RMS residual and solve tolerance
      real(wp), intent(in) :: rmse, tolerance
      !> Whether the closure changed
      logical, intent(out) :: switched

      switched = .not. self%at_target() .and. rmse < max(self%rmse, tolerance)
      if (switched) call select_rung(self, closure, self%rung + 1)
   end subroutine ramp_step

   !> Whether the current rung is the requested closure
   !>
   !> @param[in] self ramp
   pure logical function ramp_at_target(self) result(at_target)
      !> Ramp
      class(moist_model_moz_closure_ramp_type), intent(in) :: self

      at_target = self%rung == size(self%rungs)
   end function ramp_at_target

   !> Name of the current closure
   !>
   !> @param[in] self ramp
   pure function ramp_name(self) result(name)
      !> Ramp
      class(moist_model_moz_closure_ramp_type), intent(in) :: self
      !> Closure name
      character(len=:), allocatable :: name

      name = trim(self%names(self%rung))
   end function ramp_name

   !> Table note on the last step, e.g. "Closure: KH -> PSE2 at RMS residual 9.81e-06"
   !>
   !> @param[in] self ramp above its first rung
   !> @param[in] rmse RMS residual that triggered the step
   function ramp_note(self, rmse) result(text)
      !> Ramp
      class(moist_model_moz_closure_ramp_type), intent(in) :: self
      !> RMS residual
      real(wp), intent(in) :: rmse
      !> Note text
      character(len=:), allocatable :: text
      character(len=16) :: buf

      write (buf, "(es9.2)") rmse
      text = "Closure: "//trim(self%names(max(1, self%rung - 1)))//" -> "//trim(self%names(self%rung))// &
         & " at RMS residual "//trim(adjustl(buf))
   end function ramp_note

   !> Copy the closure of a rung into the closure of the solve
   !>
   !> @param[inout] self    ramp
   !> @param[inout] closure closure of the solve
   !> @param[in]    rung    rung to select
   subroutine select_rung(self, closure, rung)
      !> Ramp
      type(moist_model_moz_closure_ramp_type), intent(inout) :: self
      !> Closure of the solve
      class(moist_model_moz_closure_type), allocatable, intent(inout) :: closure
      !> Rung to select
      integer, intent(in) :: rung

      if (rung == self%rung .and. allocated(closure)) return
      if (allocated(closure)) deallocate (closure)
      allocate (closure, source=self%rungs(rung)%closure)
      self%rung = rung
   end subroutine select_rung

end module moist_model_moz_closure_ramp
