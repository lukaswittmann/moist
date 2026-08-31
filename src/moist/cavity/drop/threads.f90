!> Shared OpenMP framework for the DROP parallelization
!>
!> Routines where individual workers work independently on parts of the grid
!> require two things:
!>
!> 1. `drop_worker_slots_type`: per-thread evaluator (clone)
!>    The LSF evaluator caches screened derivatives inside itself; each
!>    thread gets its own `lsf_model` (with given highest needed derivative);
!>    if projection objective is needed, they also get a per-thread `phi`
!>
!> 2. `drop_abort_latch_type`: a shared error abort mechanism
!>    Failures insode `!$omp do` cannot be returned as easily directly;
!>    Save it inside the latch and drain the loop; after the loop the failure
!>    at the lowest grid point is raised, so the report is schedule-independent
!>
!> Usage:
!>   ```fortran
!>   call slots%init(ctx, self%lsf_model, max_deriv)
!>   call latch%reset()
!>   !$omp parallel num_threads(slots%nthreads) default(shared) ...
!>   ```
module moist_cavity_drop_threads
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_cavity_drop_objective_phi, only: moist_cavity_drop_objective_phi_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type

   implicit none (type, external)
   private

   public :: drop_worker_slots_type, drop_abort_latch_type

   !> Per-thread evaluator clones for one parallel grid loop
   !>
   !> `nthreads` is the single source of truth for the team size: it is read
   !> from the shared context and **must** be pinned on the parallel region
   !> with `num_threads(slots%nthreads)`, or a larger live team indexes the
   !> slot arrays out of bounds.
   type :: drop_worker_slots_type
      !> Per-thread LSF evaluators, one per slot
      type(lsf_thread_slot), allocatable :: lsf(:)
      !> Per-thread projection objectives; unallocated unless requested
      type(moist_cavity_drop_objective_phi_type), allocatable :: phi(:)
      !> Team size the slots were built for
      integer :: nthreads = 0
   contains
      !> Clone the evaluators for a parallel region
      procedure :: init => worker_slots_init
   end type drop_worker_slots_type

   !> Lowest-grid-point-wins error latch for a parallel grid loop
   !>
   !> A latch is written only on a failure path, so the serialization it costs
   !> is irrelevant.
   type :: drop_abort_latch_type
      !> Whether any thread has latched a failure
      logical :: requested = .false.
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
   end type drop_abort_latch_type

contains

!* ================================================================================= *!
!*                                   Worker slots                                    *!
!* ================================================================================= *!

   !> Build per-thread evaluator clones for a parallel grid loop
   !>
   !> `phi` is allocated only when `param`, `mol` and `radii` are all supplied;
   !> the projection and curvature loops do not evaluate the objective and pass
   !> none of them.
   !>
   !> @param[out] self      Worker slots
   !> @param[in]  ctx       Shared context; sole source of the team size
   !> @param[in]  lsf_model LSF template to clone
   !> @param[in]  max_deriv Highest derivative order the clones must serve
   !> @param[in]  param     DROP parameters, for the per-thread objectives
   !> @param[in]  mol       Molecular structure, for the per-thread objectives
   !> @param[in]  radii     Atomic radii, for the per-thread objectives
   subroutine worker_slots_init(self, ctx, lsf_model, max_deriv, param, mol, radii)
      !> Worker slots
      class(drop_worker_slots_type), intent(out) :: self
      !> Shared context
      type(moist_context_type), intent(in) :: ctx
      !> LSF template
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_model
      !> Highest derivative order required
      integer, intent(in) :: max_deriv
      !> DROP parameters
      type(moist_cavity_drop_parameters_type), intent(in), optional :: param
      !> Molecular structure
      type(structure_type), intent(in), optional :: mol
      !> Atomic radii
      real(wp), intent(in), optional :: radii(:)

      !> Slot index
      integer :: islot
      !> Whether the per-thread objectives are requested
      logical :: want_phi

      want_phi = present(param) .and. present(mol) .and. present(radii)

      self%nthreads = ctx%get_num_threads()
      allocate (self%lsf(self%nthreads))
      if (want_phi) allocate (self%phi(self%nthreads))

      do islot = 1, self%nthreads
         allocate (self%lsf(islot)%lsf, source=lsf_model)
         ! The screened-derivative cache is sized here, before the first
         ! %prepare call; a later upgrade would have to reallocate it.
         call self%lsf(islot)%lsf%set_max_deriv(max_deriv)
         if (want_phi) then
            call self%phi(islot)%set_parameters(param)
            call self%phi(islot)%set_input(mol, radii)
         end if
      end do
   end subroutine worker_slots_init

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
   !> worker can simply drop out of the loop afterwards.
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
         self%requested = .true.
         self%igrid = candidate
         self%status = 0
         ! `fatal_error` is intent(out) in its error argument, so a message
         ! latched over an earlier one releases it.
         call fatal_error(self%error, message)
      end if
      !$omp end critical (drop_abort_latch)
   end subroutine abort_latch_message

   !> Latch a failure described by a status code
   !>
   !> The code is stored verbatim; turning it into a diagnostic needs the
   !> caller's context and happens after the parallel region.
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
         self%requested = .true.
         self%igrid = igrid
         self%status = status
         ! A status latched over an error object: the two representations are
         ! alternatives, so the loser must go.
         if (allocated(self%error)) deallocate (self%error)
      end if
      !$omp end critical (drop_abort_latch)
   end subroutine abort_latch_status

   !> Whether a candidate failure displaces the one currently latched
   !>
   !> Call only from inside the `drop_abort_latch` critical region.
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
