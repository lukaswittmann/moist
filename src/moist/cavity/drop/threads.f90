!> Shared OpenMP framework for the DROP parallelization
!>
!> Routines where individual workers work independently on parts of the grid
!> require two things:
!>
!> 1. `drop_worker_slots_type`: per-thread evaluator (clone)
!>    The LSF evaluator caches screened derivatives inside itself; each
!>    thread gets its own worker of `lsf_model` (with given highest needed
!>    derivative); if projection objective is needed, they also get a
!>    per-thread `phi`
!>    - sized serially by `init`, filled and released by each thread for its
!>      own slot inside the parallel region, so the clones are made in
!>      parallel and their memory is first touched by the thread using it
!>    - a worker may borrow read-only data of `lsf_model`
!>      ([[moist_cavity_drop_lsf_type:clone_worker]]); the model must stay
!>      allocated and unchanged until the slot is released
!>
!> 2. `drop_abort_latch_type`: a shared error abort mechanism
!>    Failures inside `!$omp do` cannot be returned as easily directly;
!>    Save it inside the latch and drain the loop; after the loop the failure
!>    at the lowest grid point is raised, so the report is schedule-independent
!>
!> Usage:
!>   ```fortran
!>   call slots%init(ctx)
!>   call latch%reset()
!>   !$omp parallel num_threads(slots%nthreads) default(shared) private(islot, ...)
!>   islot = omp_get_thread_num() + 1
!>   call slots%fill(islot, self%lsf_model, max_deriv)
!>   ...
!>   call slots%release(islot)
!>   !$omp end parallel
!>   ```
module moist_cavity_drop_threads
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_cavity_drop_phi, only: moist_cavity_drop_phi_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type

   implicit none(type, external)
   private

   public :: drop_worker_slots_type, drop_abort_latch_type

   !> Per-thread evaluator clones for one parallel grid loop
   !>
   !> `nthreads` is the single source of truth for the team size: it is read
   !> from the shared context and **must** be pinned on the parallel region
   !> with `num_threads(slots%nthreads)`, or a larger live team indexes the
   !> slot arrays out of bounds
   type :: drop_worker_slots_type
      !> Per-thread LSF evaluators, one per slot
      type(lsf_thread_slot), allocatable :: lsf(:)
      !> Per-thread projection objectives; unallocated unless requested
      type(moist_cavity_drop_phi_type), allocatable :: phi(:)
      !> Team size the slots were built for
      integer :: nthreads = 0
   contains
      !> Size the slots for a parallel region, before entering it
      procedure :: init => worker_slots_init
      !> Clone the evaluators of one slot, from the thread owning it
      procedure :: fill => worker_slots_fill
      !> Free the evaluators of one slot, from the thread owning it
      procedure :: release => worker_slots_release
   end type drop_worker_slots_type

   !> Lowest-grid-point-wins error latch for a parallel grid loop
   !>
   !> A latch is written only on a failure path, so the serialization it costs
   !> is irrelevant
   type :: drop_abort_latch_type
      !> Whether any thread has latched a failure; read through `is_set` or `skips`
      logical, private :: requested = .false.
      !> Latched error object, when the failure carried one
      type(error_type), allocatable :: error
      !> Latched status code, when the failure carried one instead
      integer :: status = 0
      !> Grid point of the latched failure; 0 for a failure with no grid point
      integer :: igrid = 0
   contains
      !> Clear the latch before a parallel region
      procedure :: reset => abort_latch_reset
      !> Latch an existing error object, transferring ownership
      procedure :: latch_error => abort_latch_error
      !> Latch a message, building the error object here
      procedure :: latch_message => abort_latch_message
      !> Latch a status code and the grid point it came from
      procedure :: latch_status => abort_latch_status
      !> Whether any failure is latched, safe inside a parallel region
      procedure :: is_set => abort_latch_is_set
      !> Whether a loop index lies behind the latched failure
      procedure :: skips => abort_latch_skips
   end type drop_abort_latch_type

contains

!* ================================================================================= *!
!*                                   Worker slots                                    *!
!* ================================================================================= *!

   !> Size the per-thread slots of a parallel grid loop
   !>
   !> - serial and cheap: no evaluator is cloned here, each thread fills its own
   !>   slot inside the region ([[worker_slots_fill]])
   !>
   !> @param[out] self      Worker slots
   !> @param[in]  ctx       Shared context; sole source of the team size
   !> @param[in]  objective Whether the slots carry per-thread objectives, default no
   subroutine worker_slots_init(self, ctx, objective)
      !> Worker slots
      class(drop_worker_slots_type), intent(out) :: self
      !> Shared context
      type(moist_context_type), intent(in) :: ctx
      !> Whether per-thread objectives are needed
      logical, intent(in), optional :: objective

      self%nthreads = ctx%get_num_threads()
      allocate (self%lsf(self%nthreads))
      if (present(objective)) then
         if (objective) allocate (self%phi(self%nthreads))
      end if
   end subroutine worker_slots_init

   !> Clone the evaluators of one slot
   !>
   !> - called by the thread owning `islot`, inside the parallel region; slots
   !>   are disjoint, so the threads fill them concurrently
   !> - `phi` is set only when the slots carry objectives and `param`, `mol`
   !>   and `radii` are all supplied; slots made with `objective` must get all
   !>   three, or `phi` stays default and unusable
   !>
   !> @param[inout] self      Worker slots, sized by [[worker_slots_init]]
   !> @param[in]    islot     Slot of the calling thread, 1..nthreads
   !> @param[in]    lsf_model LSF template; read-only data may be borrowed
   !> @param[in]    max_deriv Highest derivative order the clone must serve
   !> @param[in]    param     DROP parameters, for the per-thread objective
   !> @param[in]    mol       Molecular structure, for the per-thread objective
   !> @param[in]    radii     Atomic radii, for the per-thread objective
   subroutine worker_slots_fill(self, islot, lsf_model, max_deriv, param, mol, radii)
      !> Worker slots
      class(drop_worker_slots_type), intent(inout) :: self
      !> Slot of the calling thread
      integer, intent(in) :: islot
      !> LSF template
      class(moist_cavity_drop_lsf_type), intent(in), target :: lsf_model
      !> Highest derivative order required
      integer, intent(in) :: max_deriv
      !> DROP parameters
      type(moist_cavity_drop_parameters_type), intent(in), optional :: param
      !> Molecular structure
      type(structure_type), intent(in), optional :: mol
      !> Atomic radii
      real(wp), intent(in), optional :: radii(:)

      call lsf_model%clone_worker(self%lsf(islot)%lsf)
      call self%lsf(islot)%lsf%set_max_deriv(max_deriv)
      if (allocated(self%phi) .and. present(param) .and. present(mol) .and. present(radii)) then
         call self%phi(islot)%set_parameters(param)
         call self%phi(islot)%set_input(mol, radii)
      end if
   end subroutine worker_slots_fill

   !> Free the evaluators of one slot
   !>
   !> - called by the thread owning `islot` before leaving the region, so the
   !>   clones are freed in parallel and not by the master after it
   !>
   !> @param[inout] self  Worker slots
   !> @param[in]    islot Slot of the calling thread
   subroutine worker_slots_release(self, islot)
      !> Worker slots
      class(drop_worker_slots_type), intent(inout) :: self
      !> Slot of the calling thread
      integer, intent(in) :: islot

      if (allocated(self%lsf(islot)%lsf)) deallocate (self%lsf(islot)%lsf)
   end subroutine worker_slots_release

!* ================================================================================= *!
!*                                    Abort latch                                    *!
!* ================================================================================= *!

   !> Clear the latch before entering a parallel region
   !>
   !> @param[inout] self Abort latch
   subroutine abort_latch_reset(self)
      !> Abort latch
      class(drop_abort_latch_type), intent(inout) :: self

      self%requested = .false.
      self%status = 0
      self%igrid = 0
      if (allocated(self%error)) deallocate (self%error)
   end subroutine abort_latch_reset

   !> Latch an error object raised by a worker, transferring ownership
   !>
   !> `err` is deallocated on return whether or not this caller won, so the
   !> worker can simply drop out of the loop afterwards
   !>
   !> @param[inout] self  Abort latch
   !> @param[inout] err   Error to latch; deallocated on return
   !> @param[in]    igrid Grid point that failed; absent means a setup failure
   subroutine abort_latch_error(self, err, igrid)
      !> Abort latch
      class(drop_abort_latch_type), intent(inout) :: self
      !> Error to latch
      type(error_type), allocatable, intent(inout) :: err
      !> Failing grid point
      integer, intent(in), optional :: igrid

      !> Grid point this candidate carries
      integer :: candidate

      candidate = 0
      if (present(igrid)) candidate = igrid

      !$omp critical (drop_abort_latch)
      if (latch_wins(self, candidate)) then
         !$omp atomic write
         self%requested = .true.
         self%igrid = candidate
         self%status = 0
         call move_alloc(err, self%error)
      end if
      !$omp end critical (drop_abort_latch)

      if (allocated(err)) deallocate (err)
   end subroutine abort_latch_error

   !> Latch a failure described by a message
   !>
   !> @param[inout] self    Abort latch
   !> @param[in]    message Diagnostic to report
   !> @param[in]    igrid   Grid point that failed; absent means a setup failure
   subroutine abort_latch_message(self, message, igrid)
      !> Abort latch
      class(drop_abort_latch_type), intent(inout) :: self
      !> Diagnostic message
      character(len=*), intent(in) :: message
      !> Failing grid point
      integer, intent(in), optional :: igrid

      !> Grid point this candidate carries
      integer :: candidate

      candidate = 0
      if (present(igrid)) candidate = igrid

      !$omp critical (drop_abort_latch)
      if (latch_wins(self, candidate)) then
         !$omp atomic write
         self%requested = .true.
         self%igrid = candidate
         self%status = 0
         ! `fatal_error` is intent(out) in its error argument, so a message
         ! latched over an earlier one releases it
         call fatal_error(self%error, message)
      end if
      !$omp end critical (drop_abort_latch)
   end subroutine abort_latch_message

   !> Latch a failure described by a status code
   !>
   !> The code is stored verbatim; turning it into a diagnostic needs the
   !> caller's context and happens after the parallel region
   !>
   !> @param[inout] self   Abort latch
   !> @param[in]    status Status code to latch
   !> @param[in]    igrid  Grid point that failed
   subroutine abort_latch_status(self, status, igrid)
      !> Abort latch
      class(drop_abort_latch_type), intent(inout) :: self
      !> Status code
      integer, intent(in) :: status
      !> Failing grid point
      integer, intent(in) :: igrid

      !$omp critical (drop_abort_latch)
      if (latch_wins(self, igrid)) then
         !$omp atomic write
         self%requested = .true.
         self%igrid = igrid
         self%status = status
         ! A status latched over an error object: the two representations are
         ! alternatives, so the loser must go
         if (allocated(self%error)) deallocate (self%error)
      end if
      !$omp end critical (drop_abort_latch)
   end subroutine abort_latch_status

   !> Whether any thread has latched a failure
   !>
   !> - atomic read against the atomic write of the latching routines, so a
   !>   worker may ask from inside the parallel loop
   !>
   !> @param[in] self  abort latch
   function abort_latch_is_set(self) result(set)
      !> Abort latch
      class(drop_abort_latch_type), intent(in) :: self
      !> Whether a failure is latched
      logical :: set

      !$omp atomic read
      set = self%requested
   end function abort_latch_is_set

   !> Whether a loop index lies behind the latched failure
   !>
   !> - `.false.` while nothing is latched, at the cost of one atomic read
   !> - the latched index is read under the latch's lock, so it belongs to the
   !>   failure seen; an index at or ahead of it is not skipped
   !> - a failure latched without an index (0) skips every index from 1 on
   !>
   !> @param[in] self  abort latch
   !> @param[in] idx   loop index, in the numbering the failures are latched with
   function abort_latch_skips(self, idx) result(skip)
      !> Abort latch
      class(drop_abort_latch_type), intent(in) :: self
      !> Loop index
      integer, intent(in) :: idx
      !> Whether the index lies behind the latched failure
      logical :: skip

      !$omp atomic read
      skip = self%requested
      if (.not. skip) return

      !$omp critical (drop_abort_latch)
      skip = idx > self%igrid
      !$omp end critical (drop_abort_latch)
   end function abort_latch_skips

   !> Whether a candidate failure displaces the one currently latched
   !>
   !> Call only from inside the `drop_abort_latch` critical region
   !>
   !> @param[in] self      Abort latch
   !> @param[in] candidate Grid point of the candidate failure
   !> @return    wins      `.true.` when the candidate must be stored
   pure function latch_wins(self, candidate) result(wins)
      !> Abort latch
      class(drop_abort_latch_type), intent(in) :: self
      !> Grid point of the candidate failure
      integer, intent(in) :: candidate
      !> Whether the candidate displaces the latched failure
      logical :: wins

      wins = .not. self%requested .or. candidate < self%igrid
   end function latch_wins

end module moist_cavity_drop_threads
