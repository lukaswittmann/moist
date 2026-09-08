!> Shared OpenMP framework for the DROP parallelization
!>
!> Routines where individual workers work independently on parts of the grid
!> require three things:
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
!> 3. `drop_point_scratch_type`: the per-thread state of one grid point
!>    One object, so the `private(...)` list carries one name instead of the
!>    two to seven dozen the traversals used to spell out by hand -- where
!>    every omitted name is a silent data race
!>
!> Usage:
!>   ```fortran
!>   call slots%init(ctx, self%lsf_model, max_deriv)
!>   call latch%reset()
!>   !$omp parallel num_threads(slots%nthreads) default(shared) &
!>   !$omp& private(thread_slot, igrid, pt, ...)
!>   call pt%init(nsph, iswig, ...)
!>   ...
!>   call pt%destroy()
!>   !$omp end parallel
!>   ```
module moist_cavity_drop_threads
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_cavity_drop_objective_phi, only: moist_cavity_drop_objective_phi_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_gaussian, only: iswig_workspace_type, moist_cavity_drop_iswig
   use moist_cavity_drop_derivatives_kernel, only: seed_status_message, drop_seed_state_type
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type, drop_n_seed_columns

   implicit none (type, external)
   private

   public :: drop_worker_slots_type, drop_abort_latch_type, drop_point_scratch_type

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
      !> Render the latched failure as the caller's error object
      procedure :: raise => abort_latch_raise
   end type drop_abort_latch_type

   !> Per-thread state of one grid point
   !>
   !> Every DROP derivative traversal opens a grid point the same way, and used
   !> to carry the result of that opening in a hand-maintained `private(...)`
   !> list of two to seven dozen names. An omitted name is a silent data race
   !> that no result check reliably catches, so the set is bundled here: a
   !> traversal declares one `type(drop_point_scratch_type) :: pt`, names it
   !> once, and inherits the whole set.
   !>
   !> Two kinds of member live here. The first is what [[drop_point_prologue]]
   !> writes at every point -- the anchor, its owner sphere, the multiplier, the
   !> level-set and objective jets, the derived seed state, the bordered
   !> factorization and the solved seed batch -- which is why the prologue takes
   !> this object instead of a dozen output arguments. The second is the buffers
   !> that have to be allocated once per thread *outside* the grid loop: an
   !> automatic array or a per-point `allocate` would move a heap operation into
   !> the hot loop.
   !>
   !> `kkt_rhs`, `lsf4_rrrr` and `vjp_pt` are optional members. They are
   !> requested at [[point_scratch_init]] and detected afterwards by their
   !> allocation status, so the prologue fills the fourth-derivative buffer and
   !> solves the standard seed batch exactly for the callers that asked.
   !>
   !> The type is concrete rather than polymorphic on purpose: gfortran 14 in
   !> this project cannot reliably copy a `class(...)` allocatable, and an
   !> OpenMP `private` copy is exactly that operation.
   type :: drop_point_scratch_type
      !> Anchor of the grid point
      real(wp) :: anchor(3) = 0.0_wp
      !> Owner sphere of the anchor
      integer :: owner_idx = 0
      !> Lagrange multiplier of the projection
      real(wp) :: lambda_val = 0.0_wp
      !> Level-set gradient and Hessian at the projected point
      real(wp) :: lsf1_r(3) = 0.0_wp, lsf2_rr(3, 3) = 0.0_wp
      !> Objective gradient at the projected point
      real(wp) :: phi1_r(3) = 0.0_wp
      !> Per-grid point sensitivity kernel state
      type(drop_seed_state_type) :: state
      !> Factorization reused by every solve at this grid point
      type(drop_kkt_factor_type) :: kkt_fac
      !> Solved standard batch of jet and anchor seeds; requested member
      real(wp), allocatable :: kkt_rhs(:, :)
      !> Third-derivative buffer of the level set
      real(wp), allocatable :: lsf3_rrr(:, :, :)
      !> Fourth-derivative buffer of the level set; requested member
      real(wp), allocatable :: lsf4_rrrr(:, :, :, :)
      !> Active atoms of the level set at this point, and how many there are
      integer :: n_active = 0
      integer, allocatable :: active_idx(:)
      !> Jet-contracted nuclear partials, one column per active atom; requested
      real(wp), allocatable :: vjp_pt(:, :)
      !> iSwiG neighbour cache and the sparse switching rows it feeds
      type(iswig_workspace_type) :: iswig_work
      real(wp), allocatable :: swi_rows(:, :)
   contains
      !> Allocate the per-thread buffers, on the thread that will own them
      procedure :: init => point_scratch_init
      !> Release them again
      procedure :: destroy => point_scratch_destroy
   end type drop_point_scratch_type

contains

!* ================================================================================= *!
!*                                   Worker slots                                    *!
!* ================================================================================= *!

   !> Build per-thread evaluator clones for a parallel grid loop
   !>
   !> `phi` is allocated only when `param`, `mol` and `radii` are all supplied;
   !> the projection and curvature loops do not evaluate the objective and pass
   !> none of them
   !>
   !> @param[out] self      Worker slots
   !> @param[in]  ctx       Shared context; sole source of the team size
   !> @param[in]  lsf_model LSF template to clone
   !> @param[in]  max_deriv Highest derivative order the clones must serve
   !> @param[in]  param     DROP parameters, for the per-thread objectives
   !> @param[in]  mol       Molecular structure, for the per-thread objectives
   !> @param[in]  radii     Atomic radii, for the per-thread objectives
   !> @param[in]  nthreads  Team size to build for, overriding the context's
   subroutine worker_slots_init(self, ctx, lsf_model, max_deriv, param, mol, radii, nthreads)
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
      !> Team size to build for; the context's own is used when absent
      integer, intent(in), optional :: nthreads

      !> Slot index
      integer :: islot
      !> Whether the per-thread objectives are requested
      logical :: want_phi

      want_phi = present(param) .and. present(mol) .and. present(radii)

      ! A serial caller passes `nthreads = 1` rather than cloning the LSF -- and
      ! with it the screened-derivative cache -- once per thread of a team it
      ! never starts. Every parallel caller leaves it absent, so the context
      ! stays the single source of truth for a team that is actually run.
      self%nthreads = ctx%get_num_threads()
      if (present(nthreads)) self%nthreads = nthreads
      allocate (self%lsf(self%nthreads))
      if (want_phi) allocate (self%phi(self%nthreads))

      do islot = 1, self%nthreads
         allocate (self%lsf(islot)%lsf, source=lsf_model)
         ! The screened-derivative cache is sized here, before the first
         ! %prepare call; a later upgrade would have to reallocate it
         call self%lsf(islot)%lsf%set_max_deriv(max_deriv)
         if (want_phi) then
            call self%phi(islot)%set_parameters(param)
            call self%phi(islot)%set_input(mol, radii)
         end if
      end do
   end subroutine worker_slots_init

!* ================================================================================= *!
!*                                   Point scratch                                   *!
!* ================================================================================= *!

   !> Allocate one thread's per-point buffers
   !>
   !> Call from inside the parallel region, on the thread that will own the
   !> scratch, so the storage is first touched by the thread that writes it --
   !> the same rule the worker slots and the sparse Hessian accumulators follow.
   !>
   !> The three optional members are allocated only when asked for. Nothing else
   !> distinguishes them: a caller that does not request `kkt_rhs` gets a
   !> prologue that does not solve the standard seed batch, and one that does
   !> not request `lsf4_rrrr` gets a prologue that never asks the level set for
   !> its fourth derivative.
   !>
   !> @param[out] self            Point scratch of one thread
   !> @param[in]  nsph            Atoms the per-atom buffers are sized for
   !> @param[in]  iswig           Switching function whose geometry is cached
   !> @param[in]  want_seed_batch Whether the standard seed batch is solved
   !> @param[in]  want_lsf4       Whether the fourth level-set derivative is read
   !> @param[in]  want_vjp        Whether a jet-contracted nuclear row is needed
   subroutine point_scratch_init(self, nsph, iswig, want_seed_batch, want_lsf4, want_vjp)
      !> Point scratch of one thread
      class(drop_point_scratch_type), intent(out) :: self
      !> Atoms the per-atom buffers are sized for
      integer, intent(in) :: nsph
      !> Switching function whose geometry the neighbour cache will hold
      class(moist_cavity_drop_iswig), intent(in) :: iswig
      !> Whether the standard seed batch is solved at every point
      logical, intent(in), optional :: want_seed_batch
      !> Whether the fourth level-set derivative is read at every point
      logical, intent(in), optional :: want_lsf4
      !> Whether a jet-contracted nuclear row buffer is needed
      logical, intent(in), optional :: want_vjp

      allocate (self%lsf3_rrr(3, 3, 3), source=0.0_wp)
      allocate (self%active_idx(nsph))
      call self%iswig_work%init(iswig)
      ! Sized to `nsph` rather than to the workspace capacity: `n_nb` is bounded
      ! by the atom count on every traversal, so this stays valid even if the
      ! workspace has to grow itself.
      allocate (self%swi_rows(3, nsph))

      if (present(want_seed_batch)) then
         if (want_seed_batch) then
            allocate (self%kkt_rhs(4, drop_n_seed_columns), source=0.0_wp)
         end if
      end if
      if (present(want_lsf4)) then
         if (want_lsf4) allocate (self%lsf4_rrrr(3, 3, 3, 3), source=0.0_wp)
      end if
      if (present(want_vjp)) then
         if (want_vjp) allocate (self%vjp_pt(3, nsph), source=0.0_wp)
      end if
   end subroutine point_scratch_init

   !> Release one thread's per-point buffers
   !>
   !> @param[inout] self Point scratch of one thread
   subroutine point_scratch_destroy(self)
      !> Point scratch of one thread
      class(drop_point_scratch_type), intent(inout) :: self

      if (allocated(self%kkt_rhs)) deallocate (self%kkt_rhs)
      if (allocated(self%lsf3_rrr)) deallocate (self%lsf3_rrr)
      if (allocated(self%lsf4_rrrr)) deallocate (self%lsf4_rrrr)
      if (allocated(self%active_idx)) deallocate (self%active_idx)
      if (allocated(self%vjp_pt)) deallocate (self%vjp_pt)
      if (allocated(self%swi_rows)) deallocate (self%swi_rows)
      call self%iswig_work%destroy()
      self%n_active = 0
   end subroutine point_scratch_destroy

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
         self%requested = .true.
         self%igrid = igrid
         self%status = status
         ! A status latched over an error object: the two representations are
         ! alternatives, so the loser must go
         if (allocated(self%error)) deallocate (self%error)
      end if
      !$omp end critical (drop_abort_latch)
   end subroutine abort_latch_status

   !> Render the latched failure as an error object for the caller
   !>
   !> The epilogue every parallel grid traversal shares. A failure reaches the
   !> latch in one of two shapes and leaves in one: an LSF or KKT failure
   !> arrives as a ready-made error object and is handed straight on, while a
   !> kernel degeneracy arrives as a bare status code and needs the traversal's
   !> own name -- which the latch never sees -- to become a diagnostic.
   !>
   !> Nothing is raised when no failure was latched, so the caller's `error`
   !> stays unallocated and a `call raise` on a clean latch is a no-op. The
   !> caller keeps whatever unwinding of its own the failure needs -- stopping
   !> its timer, restoring an accumulator it promised to leave untouched -- and
   !> returns; this routine only builds the error.
   !>
   !> @param[inout] self    Abort latch; a latched error object is moved out
   !> @param[in]    context Calling routine, used to prefix the diagnostic
   !> @param[out]   error   Error object, allocated when a failure was latched
   subroutine abort_latch_raise(self, context, error)
      !> Abort latch
      class(drop_abort_latch_type), intent(inout) :: self
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Rendered grid index of a degenerate point; fixed length, so the message
      !> build is thread safe
      character(len=32) :: idx

      if (.not. self%requested) return

      if (allocated(self%error)) then
         call move_alloc(self%error, error)
      else
         write (idx, "(i0)") self%igrid
         call fatal_error(error, context//": "//seed_status_message(self%status)// &
                          " at grid point "//trim(idx))
      end if
   end subroutine abort_latch_raise

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
