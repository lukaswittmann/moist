!> Structure-of-arrays batch of projection points for the batched DROP solvers
!>
!> - one column per point: iterate, multiplier, level set jet, objective, status
!> - one run per anchor: a contiguous column range plus the anchor record
!> - flush-when-full accumulator: [[batch_would_overflow]] asks before an anchor's
!>   points are appended, so a run is never split across flushes, and a run longer
!>   than the target is accepted as one oversize batch
module moist_cavity_drop_projection_batch
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   implicit none(type, external)
   private

   public :: projection_batch_type
   public :: point_pending, point_converged, point_accepted, point_captured
   public :: point_failed_iterations, point_failed_singular, point_failed_gradient
   public :: point_failed_not_kkt, point_failed_saddle, point_failed_unreached
   public :: point_status_label, default_batch_size

   !* ================================================================================= *!
   !*                                Column status codes                                *!
   !* ================================================================================= *!

   ! TODO: maybe there is a more elegant way to handle this

   !> Not yet solved, or reopened for another Newton pass
   integer, parameter :: point_pending = 0
   !> Newton converged; optimality not yet verified
   integer, parameter :: point_converged = 1
   !> Converged and verified as a constrained local minimum
   integer, parameter :: point_accepted = 2
   !> Duplicate: Newton prediction on an earlier accepted column of its run
   integer, parameter :: point_captured = 3
   !> Newton ran out of iterations
   integer, parameter :: point_failed_iterations = -1
   !> Singular KKT system
   integer, parameter :: point_failed_singular = -2
   !> Vanishing level set gradient
   integer, parameter :: point_failed_gradient = -4
   !> Converged point violates the KKT conditions
   integer, parameter :: point_failed_not_kkt = -5
   !> Saddle the Riemannian escape could not leave
   integer, parameter :: point_failed_saddle = -6
   !> Iterate beyond the level set's screening reach
   integer, parameter :: point_failed_unreached = -7

   !> Default flush target (columns per batch)
   integer, parameter :: default_batch_size = 256

   !> Batch of projection points, structure of arrays
   !>
   !> Column arrays are sized `capacity`, valid on `1:n`; run arrays are sized
   !> `run_capacity`, valid on `1:nrun`
   type :: projection_batch_type
      !> Flush target; a batch holding points flushes before exceeding it
      integer :: target_size = default_batch_size
      !> Columns in use
      integer :: n = 0
      !> Allocated columns
      integer :: capacity = 0

      !> Iterate (3, capacity)
      real(wp), allocatable :: x(:, :)
      !> Multiplier
      real(wp), allocatable :: lambda(:)
      !> Level set value at `x`
      real(wp), allocatable :: S(:)
      !> Level set gradient at `x` (3, capacity)
      real(wp), allocatable :: gradS(:, :)
      !> Level set Hessian at `x`, packed (6, capacity) as the LSF returns it
      real(wp), allocatable :: hessS(:, :)
      !> Anchor objective `phi_alpha/2 ||x - anchor||^2`
      real(wp), allocatable :: phi(:)
      !> Anchor distance
      real(wp), allocatable :: rho(:)
      !> Unit surface normal (3, capacity)
      real(wp), allocatable :: normal(:, :)
      !> Column status, one of the `point_*` codes
      integer, allocatable :: status(:)
      !> Run index of the column
      integer, allocatable :: run_of(:)

      !> Runs in use
      integer :: nrun = 0
      !> Allocated runs
      integer :: run_capacity = 0
      !> First column of the run
      integer, allocatable :: run_start(:)
      !> Columns in the run
      integer, allocatable :: run_len(:)
      !> Caller's anchor index
      integer, allocatable :: run_anchor(:)
      !> Owning atom of the anchor
      integer, allocatable :: run_owner(:)
      !> Anchor position (3, run_capacity)
      real(wp), allocatable :: run_xyz(:, :)
      !> Level set value at the anchor
      real(wp), allocatable :: run_S(:)
      !> Level set gradient at the anchor (3, run_capacity)
      real(wp), allocatable :: run_gradS(:, :)
      !> Whether `run_S` and `run_gradS` were recorded
      logical, allocatable :: run_jet(:)
   contains
      !> Set the flush target and presize the storage
      procedure :: init => batch_init
      !> Forget all columns and runs, keeping the storage
      procedure :: clear => batch_clear
      !> Whether appending `n_next` points must be preceded by a flush
      procedure :: would_overflow => batch_would_overflow
      !> Open the run of one anchor
      procedure :: begin_run => batch_begin_run
      !> Append a point to the open run
      procedure :: append => batch_append
      !> Append a copy of another batch's solved column to the open run
      procedure :: append_copy => batch_append_copy
      !> Record the level set jet at a run's anchor
      procedure :: set_anchor_jet => batch_set_anchor_jet
      !> Level set jet at a run's anchor; aborts when never recorded
      procedure :: anchor_jet => batch_anchor_jet
      !> Grow the column storage
      procedure :: reserve => batch_reserve
      !> Grow the run storage
      procedure :: reserve_runs => batch_reserve_runs
      !> Release all storage
      procedure :: destroy => batch_destroy
   end type projection_batch_type

contains

   !> Human-readable label of a column status
   !>
   !> @param[in] status Column status code
   pure function point_status_label(status) result(label)
      !> Column status code
      integer, intent(in) :: status
      !> Label
      character(len=:), allocatable :: label

      select case (status)
      case (point_pending)
         label = "pending"
      case (point_converged)
         label = "converged"
      case (point_accepted)
         label = "accepted"
      case (point_captured)
         label = "duplicate of an accepted foot"
      case (point_failed_iterations)
         label = "Newton iteration budget exhausted"
      case (point_failed_singular)
         label = "singular KKT system"
      case (point_failed_gradient)
         label = "vanishing level set gradient"
      case (point_failed_not_kkt)
         label = "not a KKT point"
      case (point_failed_saddle)
         label = "unescapable saddle"
      case (point_failed_unreached)
         label = "beyond the level set reach"
      case default
         label = "unknown"
      end select
   end function point_status_label

   !* ================================================================================= *!
   !*                                  Batch lifecycle                                  *!
   !* ================================================================================= *!

   !> Set the flush target and presize the storage
   !>
   !> @param[inout] self        Batch
   !> @param[in]    target_size Flush target (columns), positive
   !> @param[out]   error       Non-positive target
   !> @param[in]    capacity    Initial column capacity (default `target_size`)
   subroutine batch_init(self, target_size, error, capacity)
      class(projection_batch_type), intent(inout) :: self
      integer, intent(in) :: target_size
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: capacity

      integer :: cap

      if (target_size < 1) then
         call fatal_error(error, "Projection batch: target size must be positive")
         return
      end if
      self%target_size = target_size
      call self%clear()
      cap = target_size
      if (present(capacity)) cap = max(1, capacity)
      call self%reserve(cap)
      call self%reserve_runs(max(1, cap/4))
   end subroutine batch_init

   !> Forget all columns and runs, keeping the storage
   !>
   !> @param[inout] self Batch
   subroutine batch_clear(self)
      class(projection_batch_type), intent(inout) :: self

      self%n = 0
      self%nrun = 0
   end subroutine batch_clear

   !> Whether appending `n_next` points must be preceded by a flush
   !>
   !> - an empty batch never overflows, so a run longer than the target is
   !>   accepted whole rather than split or refused
   !>
   !> @param[in] self   Batch
   !> @param[in] n_next Points the next anchor will append
   pure function batch_would_overflow(self, n_next) result(overflow)
      class(projection_batch_type), intent(in) :: self
      integer, intent(in) :: n_next
      logical :: overflow

      overflow = self%n > 0 .and. self%n + n_next > self%target_size
   end function batch_would_overflow

   !> Release all storage
   !>
   !> @param[inout] self Batch
   subroutine batch_destroy(self)
      class(projection_batch_type), intent(inout) :: self

      if (allocated(self%x)) deallocate (self%x)
      if (allocated(self%lambda)) deallocate (self%lambda)
      if (allocated(self%S)) deallocate (self%S)
      if (allocated(self%gradS)) deallocate (self%gradS)
      if (allocated(self%hessS)) deallocate (self%hessS)
      if (allocated(self%phi)) deallocate (self%phi)
      if (allocated(self%rho)) deallocate (self%rho)
      if (allocated(self%normal)) deallocate (self%normal)
      if (allocated(self%status)) deallocate (self%status)
      if (allocated(self%run_of)) deallocate (self%run_of)
      if (allocated(self%run_start)) deallocate (self%run_start)
      if (allocated(self%run_len)) deallocate (self%run_len)
      if (allocated(self%run_anchor)) deallocate (self%run_anchor)
      if (allocated(self%run_owner)) deallocate (self%run_owner)
      if (allocated(self%run_xyz)) deallocate (self%run_xyz)
      if (allocated(self%run_S)) deallocate (self%run_S)
      if (allocated(self%run_gradS)) deallocate (self%run_gradS)
      if (allocated(self%run_jet)) deallocate (self%run_jet)
      self%n = 0
      self%capacity = 0
      self%nrun = 0
      self%run_capacity = 0
   end subroutine batch_destroy

   !* ================================================================================= *!
   !*                                  Runs and columns                                 *!
   !* ================================================================================= *!

   !> Open the run of one anchor; later appends land in it
   !>
   !> @param[inout] self       Batch
   !> @param[in]    anchor     Caller's anchor index
   !> @param[in]    anchor_xyz Anchor position
   !> @param[in]    owner      Owning atom
   subroutine batch_begin_run(self, anchor, anchor_xyz, owner)
      class(projection_batch_type), intent(inout) :: self
      integer, intent(in) :: anchor
      real(wp), intent(in) :: anchor_xyz(3)
      integer, intent(in) :: owner

      if (self%nrun + 1 > self%run_capacity) call self%reserve_runs(self%nrun + 1)
      self%nrun = self%nrun + 1
      self%run_start(self%nrun) = self%n + 1
      self%run_len(self%nrun) = 0
      self%run_anchor(self%nrun) = anchor
      self%run_owner(self%nrun) = owner
      self%run_xyz(:, self%nrun) = anchor_xyz
      self%run_jet(self%nrun) = .false.
   end subroutine batch_begin_run

   !> Append a point to the open run
   !>
   !> - the jet, objective and normal start as NaN, so a column read before
   !>   its first evaluation cannot pass for a computed value
   !>
   !> @param[inout] self   Batch
   !> @param[in]    x      Initial iterate
   !> @param[in]    lambda Initial multiplier
   !> @param[out]   icol   Column index of the new point
   subroutine batch_append(self, x, lambda, icol)
      class(projection_batch_type), intent(inout) :: self
      real(wp), intent(in) :: x(3)
      real(wp), intent(in) :: lambda
      integer, intent(out) :: icol

      real(wp) :: nan

      if (self%nrun < 1) error stop "moist DROP projection batch: append before begin_run"
      if (self%n + 1 > self%capacity) call self%reserve(self%n + 1)

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      self%n = self%n + 1
      icol = self%n
      self%x(:, icol) = x
      self%lambda(icol) = lambda
      self%S(icol) = nan
      self%gradS(:, icol) = nan
      self%hessS(:, icol) = nan
      self%phi(icol) = nan
      self%rho(icol) = nan
      self%normal(:, icol) = nan
      self%status(icol) = point_pending
      self%run_of(icol) = self%nrun
      self%run_len(self%nrun) = self%run_len(self%nrun) + 1
   end subroutine batch_append

   !> Append a copy of a solved column of `src` to the open run
   !>
   !> - iterate, multiplier, jet, objective, normal and status are copied, so
   !>   an accepted foot is not solved again
   !>
   !> @param[inout] self Batch with an open run
   !> @param[in]    src  Source batch
   !> @param[in]    j    Source column
   !> @param[out]   icol Column index of the copy
   subroutine batch_append_copy(self, src, j, icol)
      class(projection_batch_type), intent(inout) :: self
      type(projection_batch_type), intent(in) :: src
      integer, intent(in) :: j
      integer, intent(out) :: icol

      call self%append(src%x(:, j), src%lambda(j), icol)
      self%S(icol) = src%S(j)
      self%gradS(:, icol) = src%gradS(:, j)
      self%hessS(:, icol) = src%hessS(:, j)
      self%phi(icol) = src%phi(j)
      self%rho(icol) = src%rho(j)
      self%normal(:, icol) = src%normal(:, j)
      self%status(icol) = src%status(j)
   end subroutine batch_append_copy

   !> Record the level set jet at a run's anchor
   !>
   !> @param[inout] self Batch
   !> @param[in]    irun Run index
   !> @param[in]    S    Level set value at the anchor
   !> @param[in]    grad Level set gradient at the anchor
   subroutine batch_set_anchor_jet(self, irun, S, grad)
      class(projection_batch_type), intent(inout) :: self
      integer, intent(in) :: irun
      real(wp), intent(in) :: S
      real(wp), intent(in) :: grad(3)

      self%run_S(irun) = S
      self%run_gradS(:, irun) = grad
      self%run_jet(irun) = .true.
   end subroutine batch_set_anchor_jet

   !> Level set jet at a run's anchor
   !>
   !> - aborts when no evaluation at the anchor was recorded
   !>
   !> @param[in]  self Batch
   !> @param[in]  irun Run index
   !> @param[out] S    Level set value at the anchor
   !> @param[out] grad Level set gradient at the anchor
   subroutine batch_anchor_jet(self, irun, S, grad)
      class(projection_batch_type), intent(in) :: self
      integer, intent(in) :: irun
      real(wp), intent(out) :: S
      real(wp), intent(out) :: grad(3)

      if (.not. self%run_jet(irun)) then
         error stop "moist DROP projection batch: anchor jet read before it was evaluated"
      end if
      S = self%run_S(irun)
      grad = self%run_gradS(:, irun)
   end subroutine batch_anchor_jet

   !* ================================================================================= *!
   !*                                   Storage growth                                  *!
   !* ================================================================================= *!

   !> Geometric capacity at least `required`
   !>
   !> @param[in] current  Current capacity
   !> @param[in] required Minimum capacity
   pure function grown_capacity(current, required) result(capacity)
      integer, intent(in) :: current
      integer, intent(in) :: required
      integer :: capacity

      capacity = max(16, 2*current, required)
   end function grown_capacity

   !> Grow the column storage to at least `required`, keeping `1:n`
   !>
   !> @param[inout] self     Batch
   !> @param[in]    required Minimum column capacity
   subroutine batch_reserve(self, required)
      class(projection_batch_type), intent(inout) :: self
      integer, intent(in) :: required

      integer :: cap

      if (required <= self%capacity) return
      cap = grown_capacity(self%capacity, required)
      call grow_real_2d(self%x, 3, cap, self%n)
      call grow_real_1d(self%lambda, cap, self%n)
      call grow_real_1d(self%S, cap, self%n)
      call grow_real_2d(self%gradS, 3, cap, self%n)
      call grow_real_2d(self%hessS, 6, cap, self%n)
      call grow_real_1d(self%phi, cap, self%n)
      call grow_real_1d(self%rho, cap, self%n)
      call grow_real_2d(self%normal, 3, cap, self%n)
      call grow_int_1d(self%status, cap, self%n)
      call grow_int_1d(self%run_of, cap, self%n)
      self%capacity = cap
   end subroutine batch_reserve

   !> Grow the run storage to at least `required`, keeping `1:nrun`
   !>
   !> @param[inout] self     Batch
   !> @param[in]    required Minimum run capacity
   subroutine batch_reserve_runs(self, required)
      class(projection_batch_type), intent(inout) :: self
      integer, intent(in) :: required

      integer :: cap

      if (required <= self%run_capacity) return
      cap = grown_capacity(self%run_capacity, required)
      call grow_int_1d(self%run_start, cap, self%nrun)
      call grow_int_1d(self%run_len, cap, self%nrun)
      call grow_int_1d(self%run_anchor, cap, self%nrun)
      call grow_int_1d(self%run_owner, cap, self%nrun)
      call grow_real_2d(self%run_xyz, 3, cap, self%nrun)
      call grow_real_1d(self%run_S, cap, self%nrun)
      call grow_real_2d(self%run_gradS, 3, cap, self%nrun)
      call grow_logical_1d(self%run_jet, cap, self%nrun)
      self%run_capacity = cap
   end subroutine batch_reserve_runs

   !> Reallocate a real vector to `cap`, keeping its first `keep` entries
   !>
   !> @param[inout] a    Vector to grow
   !> @param[in]    cap  New length
   !> @param[in]    keep Leading entries to keep
   subroutine grow_real_1d(a, cap, keep)
      !> Vector to grow
      real(wp), allocatable, intent(inout) :: a(:)
      !> New length and leading entries to keep
      integer, intent(in) :: cap, keep
      real(wp), allocatable :: tmp(:)

      allocate (tmp(cap))
      if (allocated(a) .and. keep > 0) tmp(1:keep) = a(1:keep)
      call move_alloc(tmp, a)
   end subroutine grow_real_1d

   !> Reallocate a real matrix to `(lead, cap)`, keeping its first `keep` columns
   !>
   !> @param[inout] a    Matrix to grow
   !> @param[in]    lead Leading dimension
   !> @param[in]    cap  New column count
   !> @param[in]    keep Leading columns to keep
   subroutine grow_real_2d(a, lead, cap, keep)
      !> Matrix to grow
      real(wp), allocatable, intent(inout) :: a(:, :)
      !> Leading dimension, new column count and leading columns to keep
      integer, intent(in) :: lead, cap, keep
      real(wp), allocatable :: tmp(:, :)

      allocate (tmp(lead, cap))
      if (allocated(a) .and. keep > 0) tmp(:, 1:keep) = a(:, 1:keep)
      call move_alloc(tmp, a)
   end subroutine grow_real_2d

   !> Reallocate an integer vector to `cap`, keeping its first `keep` entries
   !>
   !> @param[inout] a    Vector to grow
   !> @param[in]    cap  New length
   !> @param[in]    keep Leading entries to keep
   subroutine grow_int_1d(a, cap, keep)
      !> Vector to grow
      integer, allocatable, intent(inout) :: a(:)
      !> New length and leading entries to keep
      integer, intent(in) :: cap, keep
      integer, allocatable :: tmp(:)

      allocate (tmp(cap))
      if (allocated(a) .and. keep > 0) tmp(1:keep) = a(1:keep)
      call move_alloc(tmp, a)
   end subroutine grow_int_1d

   !> Reallocate a logical vector to `cap`, keeping its first `keep` entries
   !>
   !> @param[inout] a    Vector to grow
   !> @param[in]    cap  New length
   !> @param[in]    keep Leading entries to keep
   subroutine grow_logical_1d(a, cap, keep)
      !> Vector to grow
      logical, allocatable, intent(inout) :: a(:)
      !> New length and leading entries to keep
      integer, intent(in) :: cap, keep
      logical, allocatable :: tmp(:)

      allocate (tmp(cap))
      if (allocated(a) .and. keep > 0) tmp(1:keep) = a(1:keep)
      call move_alloc(tmp, a)
   end subroutine grow_logical_1d

end module moist_cavity_drop_projection_batch
