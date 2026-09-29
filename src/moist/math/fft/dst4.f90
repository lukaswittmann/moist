!> Unnormalised DST-IV for the 1D-RISM Fourier-Bessel transform
!>
!> The transform is
!>
!>     Y_k = 2 sum_{j=0}^{n-1} X_j sin(pi (j+1/2)(k+1/2) / n),   k = 0..n-1,
!>
!> which the backend provides directly (ducc0's type-4 DST with `ortho` off is
!> exactly this convention`)
!>
!> It is its own inverse up to a factor of 2n, so one plan object serves both the
!> forward and the inverse Fourier-Bessel direction; the callers apply the differing
!> r/k weight diagonals around it
!>
!> The backend needs no plans, so `dst4_plan_type` carries only the transform
!> geometry and `dst4_work_type` is empty; neither `init` nor `new_work` can
!> fail, and a plan with an empty geometry transforms nothing
!>
!> Both types are kept because they give the radial grids their "immutable description
!> on the shared grid, scratch on the per-thread trafo" split: should a future backend
!> need real per-thread buffers again, they land here without touching a single call site
!>
module moist_math_fft_dst4
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_fft, only: moist_fft_dst4
   use, intrinsic :: iso_c_binding, only: c_int, c_double
   implicit none(type, external)
   private

   public :: dst4_plan_type, dst4_work_type

   !> Immutable description of a DST-IV over `nbatch` columns of length `npts`
   !>
   !> Carries no mutable state, so it is safe to share read-only across threads
   type :: dst4_plan_type
      !> Transform length
      integer :: npts = 0
      !> Number of columns transformed together
      integer :: nbatch = 0
   contains
      procedure :: init => dst4_plan_init
      procedure :: destroy => dst4_plan_destroy
      procedure :: new_work => dst4_plan_new_work
      generic :: execute => execute_one, execute_all
      procedure, private :: execute_one => dst4_execute_one
      procedure, private :: execute_all => dst4_execute_all
   end type dst4_plan_type

   !> Per-thread scratch for a `dst4_plan_type`
   type :: dst4_work_type
   contains
      procedure :: destroy => dst4_work_destroy
   end type dst4_work_type

contains

   !> Record the transform geometry; a non-positive length or batch count
   !> records an empty plan, whose `execute` is a no-op
   !>
   !> @param[out] self    Initialised plan
   !> @param[in]  npts    Transform length
   !> @param[in]  nbatch  Number of columns transformed together
   subroutine dst4_plan_init(self, npts, nbatch)
      !> Initialised plan
      class(dst4_plan_type), intent(out) :: self
      !> Transform length
      integer, intent(in) :: npts
      !> Number of columns transformed together
      integer, intent(in) :: nbatch

      self%npts = max(npts, 0)
      self%nbatch = max(nbatch, 0)
   end subroutine dst4_plan_init

   !> Reset the plan; idempotent
   !>
   !> @param[inout] self  Plan instance
   pure subroutine dst4_plan_destroy(self)
      !> Plan instance
      class(dst4_plan_type), intent(inout) :: self

      self%npts = 0
      self%nbatch = 0
   end subroutine dst4_plan_destroy

   !> Prepare the per-thread scratch; no-op with the current backend
   !>
   !> @param[in]  self   Plan the scratch belongs to
   !> @param[out] work   Initialised scratch
   pure subroutine dst4_plan_new_work(self, work)
      !> Plan the scratch belongs to
      class(dst4_plan_type), intent(in) :: self
      !> Initialised scratch
      type(dst4_work_type), intent(out) :: work

   end subroutine dst4_plan_new_work

   !> Release the scratch; idempotent no-op with the current backend
   !>
   !> @param[inout] self  Scratch instance
   pure subroutine dst4_work_destroy(self)
      !> Scratch instance
      class(dst4_work_type), intent(inout) :: self

   end subroutine dst4_work_destroy

   !> Apply the unnormalised DST-IV to a single column
   !>
   !> @param[in]    self  Plan instance (built with `nbatch == 1`)
   !> @param[inout] work  Per-thread scratch (unused)
   !> @param[inout] x     Input column of `npts` values; preserved
   !> @param[out]   y     Transformed column of `npts` values
   !> @param[out]   error Set if the ducc0 backend failed
   subroutine dst4_execute_one(self, work, x, y, error)
      !> Plan instance
      class(dst4_plan_type), intent(in) :: self
      !> Per-thread scratch
      type(dst4_work_type), intent(inout) :: work
      !> Input column
      real(wp), intent(inout), contiguous :: x(:)
      !> Transformed column
      real(wp), intent(out), contiguous :: y(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (self%npts < 1) return
      if (moist_fft_dst4(int(self%npts, c_int), 1_c_int, x, y, 1.0_c_double) /= 0) then
         call fatal_error(error, "dst4: DST-IV backend failed")
      end if
   end subroutine dst4_execute_one

   !> Apply the unnormalised DST-IV to every column of a batch
   !>
   !> @param[in]    self  Plan instance (built with a matching `nbatch`)
   !> @param[inout] work  Per-thread scratch (unused)
   !> @param[inout] x     Input columns, shape (npts, nbatch); preserved
   !> @param[out]   y     Transformed columns, shape (npts, nbatch)
   !> @param[out]   error Set if the ducc0 backend failed
   subroutine dst4_execute_all(self, work, x, y, error)
      !> Plan instance
      class(dst4_plan_type), intent(in) :: self
      !> Per-thread scratch
      type(dst4_work_type), intent(inout) :: work
      !> Input columns
      real(wp), intent(inout), contiguous :: x(:, :)
      !> Transformed columns
      real(wp), intent(out), contiguous :: y(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (self%npts < 1 .or. self%nbatch < 1) return
      if (moist_fft_dst4(int(self%npts, c_int), int(self%nbatch, c_int), &
         & x, y, 1.0_c_double) /= 0) then
         call fatal_error(error, "dst4: batched DST-IV backend failed")
      end if
   end subroutine dst4_execute_all

end module moist_math_fft_dst4
