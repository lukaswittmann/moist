!> Per-thread worker of the batched DROP projection
!>
!> Two passes per chunk of anchors, the seeder build between them
!> ([[projection_seeder_type]]):
!>
!> - local pass ([[worker_solve_local]], into a [[local_pass_type]]): local
!>   Newton and polish from every anchor above the switching cut; per anchor
!>   the next step (finalize, window or extended search) and the first search
!>   radius, which the build sees ([[search_geometry]])
!> - search pass ([[worker_search]]): window seeds plus the local point as one
!>   pass-2 run, flushed when full, or the extended search; batched Newton and
!>   polish, a seed landing on an accepted point of its run stopped as a
!>   duplicate; per-anchor finalize (dedupe, order, CP Jacobian, gates, branch
!>   weights, cut) into the worker's row buffer
!> - [[worker_process_chunk]]: both passes on one chunk, its own seeder built
!>   between them unless handed one (single-anchor entry, tests)
!>
!> Search radius `min(reach, window(pdist_k))`, `k` the nearest fully open point
!> found, the reach `branch%reach_off` without one:
!>
!> - beyond the reach the gate is closed; beyond the window of a point with gate
!>   one nothing has weight
!> - a partly open or closed point never shortens it
!> - level 1: a partly open point shortens it too and a grown ball restarts at
!>   `nearest_count` seeds; the search only replaces a failed local solve, so
!>   weight behind a partly open point may be missed
!>
!> Extended search (local point not fully open, or none):
!>
!> - nearest-first seeds in a ball: the window of the local point, else of the
!>   anchor's Newton length
!> - rounds batched over the chunk's anchors, seeds per anchor doubled each
!>   round, radius retaken from the fully open points accepted so far
!> - done once its seeds within the seeder's `reach` of the radius are solved
!>   and the ball covers the radius; else the ball grows (to the radius when set,
!>   else doubled up to the reach) and is seeded again: the outer shell, or
!>   without shell queries the whole ball minus solved seeds
!> - every run starts with the anchor's accepted points, so its last run is the
!>   one finalized; a round carries over only latest runs, reduced to accepted
!>   points
!>
!> Levels, by the anchor's local point:
!>
!> | level | fully open                                | partly open or closed | none (failed solve) |
!> |-------|-------------------------------------------|-----------------------|---------------------|
!> | 1     | finalized                                 | finalized             | extended search     |
!> | 2     | window; finalized at `pdist <= 0.50` Bohr | extended search       | extended search     |
!> | 3     | window                                    | extended search       | extended search     |
!> | 4     | window                                    | extended search       | extended search     |
!>
!> - zero window: finalized at every level
!> - level 2 shortcut (`level2_window_pdist`) counted in `n_window_skipped`
!> - other levels refused
!> - lattice seeder with the level's settings (`seed%init_level`): seeds
!>   pass all four filters at level 2 (back-facing one only when the horizon
!>   switch closes such normals), the tangent cone alone at level 3, none at
!>   levels 1 and 4
!> - nothing accepted after the search: `lost`, no row; also `unseeded` when no
!>   search gave a seed (no surface within the reach)
!> - `nearest_count`: starting guess, not a correctness parameter
!> - `full_reach` (debug): every search with ball and radius at the reach; same
!>   rows as the adaptive radius as far as the seeder finds the same points
!> - build-fatal, returned in `error`: LSF failure, anchor on its owner center,
!>   singular tangent block
module moist_cavity_drop_projection_worker
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type, is_valid_proj_level
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_projection_batch, only: projection_batch_type, point_accepted, point_captured
   use moist_cavity_drop_projection_parameters, only: drop_solve_parameters_type, drop_polish_parameters_type, &
                                                     drop_finalize_parameters_type
   use moist_cavity_drop_projection_newton, only: newton_scratch_type, newton_solve_batch
   use moist_cavity_drop_projection_polish, only: polish_batch, saddle_cache_type
   use moist_cavity_drop_projection_finalize, only: finalize_anchor_type, &
                                                    finalize_scratch_type, projection_rows_type, &
                                                    finalize_anchor, point_gate
   use moist_cavity_drop_projection_seed, only: projection_seeder_type, seeder_scratch_type, seed_request_type, &
                                                seed_build_type
   use moist_utils_mem, only: grow_array
   use moist_math_sorter, only: sort_index
   implicit none(type, external)
   private

   public :: projection_worker_type, local_pass_type, search_geometry, search_order

   !> Seeds per anchor in the first round of the extended search
   integer, parameter :: default_nearest_count = 4
   !> Largest `take` before doubling would overflow default integer
   integer, parameter :: max_near_take = ishft(huge(1), -1)
   !> Smallest first ball of an extended search, in units of the reach
   real(wp), parameter :: min_ball_fraction = 0.02_wp
   !> Level 2: a fully open local point at or below this distance (Bohr) is
   !> finalized without seeding its window
   real(wp), parameter :: level2_window_pdist = 0.50_wp
   !> Search-pass cost of an extended search in units of a window search, for
   !> the chunk order ([[search_order]])
   integer, parameter :: extended_cost = 16

   !> Next step of an anchor after the local pass
   !> Finalize from the local point
   integer, parameter :: next_finalize = 0
   !> Seed the branch window behind the fully open local point
   integer, parameter :: next_window = 1
   !> Extended search
   integer, parameter :: next_extended = 2

   !> Local pass of one chunk of anchors, kept for its search pass
   type :: local_pass_type
      !> Pass-1 batch: one run per anchor above the switching cut, holding its
      !> local Newton column
      type(projection_batch_type) :: batch
      !> Next step of each run, one of the `next_*` codes
      integer, allocatable :: next(:)
      !> Accepted local column of each run, 0 for none
      integer, allocatable :: jloc(:)
      !> First search radius of each run (Bohr); 0 for a run finalized from its
      !> local point
      real(wp), allocatable :: radius(:)
   contains
      !> Whether a run goes to the seeder in the search pass
      procedure :: seeded => local_seeded
      !> Whether a run's search may grow its ball beyond the first radius
      procedure :: grows => local_grows
   end type local_pass_type

   !> One anchor of the extended search
   type :: queue_entry_type
      !> Local run of the anchor
      integer :: run = 0
      !> Its accepted local column, 0 for none
      integer :: jloc = 0
      !> Its latest run in `passn`, 0 for none
      integer :: tail = 0
      !> First seed of its list in `near_seeds`
      integer :: first = 0
      !> Last seed of its list in `near_seeds`
      integer :: last = 0
      !> Next seed of its list in `near_seeds`
      integer :: next = 0
      !> Seeds of its next round
      integer :: take = 0
      !> Seeder reach of its seeds (Bohr)
      real(wp) :: reach = 0.0_wp
      !> Radius its points ask for (Bohr)
      real(wp) :: radius = 0.0_wp
      !> Radius of the ball its seeds cover (Bohr)
      real(wp) :: ball = 0.0_wp
   end type queue_entry_type

   !> Per-thread projection worker
   type :: projection_worker_type
      !> Projection level (1, 2, 3 or 4)
      integer :: level = 3
      !> Thread-local level set clone
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      !> Atom centers (3, nat)
      real(wp), allocatable :: centers(:, :)
      !> Solve-stage parameters
      type(drop_solve_parameters_type) :: solve
      !> Polish-stage parameters
      type(drop_polish_parameters_type) :: polish
      !> Finalize-stage parameters
      type(drop_finalize_parameters_type) :: finalize
      !> Solve parameters of the seeded pass: `solve` with duplicate capture
      type(drop_solve_parameters_type) :: solve_seeded
      !> Anchor objective weight, set by `init`
      real(wp) :: phi_alpha
      !> Columns per batch
      integer :: batch_size = 0
      !> Distance from which no point carries weight (Bohr): the largest search radius
      real(wp) :: reach_off = 0.0_wp
      !> Search every anchor to the reach, whatever its points (debug)
      logical :: full_reach = .false.
      !> Seeder of the level's settings, from `init`; `process_chunk` builds it
      !> for its chunk
      type(projection_seeder_type) :: seeder
      !> This thread's scratch; holds the seeds of the last request
      type(seeder_scratch_type) :: scratch
      !> Gradient-norm switch
      type(moist_math_smoothstep_sigmoid_type) :: f_crit
      !> Focal switch on the focal margin `beta_min(B)/phi_alpha`
      type(moist_math_smoothstep_sigmoid_type) :: f_foc
      !> Horizon switch on the signed `n_surf . n_sph`
      type(moist_math_smoothstep_sigmoid_type) :: f_hor
      !> Prune switch
      type(moist_math_smoothstep_sigmoid_type) :: f_wleb
      !> Reach gate on the projection distance, falling from `reach_on` to `reach_off`
      type(moist_math_smoothstep_sigmoid_type) :: f_reach
      !> Branch weight model
      type(branch_weight_type) :: branch_weight
      !> Seeds per anchor in the first round of the extended search, doubled per round
      integer :: nearest_count = default_nearest_count
      !> Extended-search batch of the current round, swapped with `passn_last` per round
      type(projection_batch_type), allocatable :: passn
      !> Extended-search batch of the last round
      type(projection_batch_type), allocatable :: passn_last
      !> Pass-2 (seeded) batch
      type(projection_batch_type) :: pass2
      !> Anchors in the extended search
      type(queue_entry_type), allocatable :: queue(:)
      !> Nearest-first seeds of those anchors (3, :)
      real(wp), allocatable :: near_seeds(:, :)
      !> Seed solved already from the anchor's smaller ball, never solved again
      logical, allocatable :: near_done(:)
      !> Newton scratch
      type(newton_scratch_type) :: nscratch
      !> Finalize scratch
      type(finalize_scratch_type) :: fscratch
      !> Persistent per-anchor saddle cache across rounds of `solve_queue`
      type(saddle_cache_type) :: ext_saddles
      !> Rows appended by every chunk this worker processed
      type(projection_rows_type) :: rows
      !> Anchors searched (window or extended search) since `init`
      integer :: n_seeded = 0
      !> Seeds solved since `init`
      integer :: n_seeds = 0
      !> Seeds of a grown ball skipped as solved from the smaller one, since `init`
      integer :: n_skipped = 0
      !> Anchors without an accepted local point, seeded since `init`
      integer :: n_escalated = 0
      !> Seeded columns stopped as duplicates of an accepted point, since `init`
      integer :: n_captured = 0
      !> Extended searches, every one whatever the anchor's local point, since `init`
      integer :: n_extended = 0
      !> Level 2: anchors finalized from a fully open local point within
      !> `level2_window_pdist`, their window not seeded, since `init`
      integer :: n_window_skipped = 0
      !> Largest search radius of an anchor since `init` (Bohr)
      real(wp) :: pdist_search_max = 0.0_wp
      !> Anchors searched to the reach from the start since `init`
      integer :: n_full = 0
   contains
      !> Configure the worker from the cavity parameters
      procedure :: init => worker_init
      !> Local pass of one chunk of anchors
      procedure :: solve_local => worker_solve_local
      !> Search pass of one chunk of anchors, from its local pass
      procedure :: search => worker_search
      !> Run both passes on one chunk of anchors
      procedure :: process_chunk => worker_process_chunk
      !> Single-anchor entry: project one anchor and return its rows
      procedure :: project_anchor => worker_project_anchor
      !> Largest distance inside the branch window behind a fully open point
      procedure :: window_radius => worker_window_radius
   end type projection_worker_type

contains

   !> Whether run `r` goes to the seeder in the search pass
   !>
   !> @param[in] self Local pass
   !> @param[in] r    Run index
   pure function local_seeded(self, r) result(seeded)
      class(local_pass_type), intent(in) :: self
      integer, intent(in) :: r
      logical :: seeded

      seeded = self%next(r) /= next_finalize
   end function local_seeded

   !> Whether the search of run `r` may grow its ball beyond the first radius,
   !> up to the reach: the extended search
   !>
   !> @param[in] self Local pass
   !> @param[in] r    Run index
   pure function local_grows(self, r) result(grows)
      class(local_pass_type), intent(in) :: self
      integer, intent(in) :: r
      logical :: grows

      grows = self%next(r) == next_extended
   end function local_grows

   !> Build input of a search pass: the geometry and the largest search
   !> ball of every anchor the pass will seed, in chunk and run order
   !>
   !> @param[in]  local    Local passes of the chunks, in chunk order
   !> @param[in]  centers  Atom centers (3, nat)
   !> @param[in]  radii    Atom radii (nat)
   !> @param[in]  reach    Largest search radius (Bohr)
   !> @param[in]  nthreads Threads of the build
   !> @param[out] geometry Build input
   subroutine search_geometry(local, centers, radii, reach, nthreads, geometry)
      type(local_pass_type), intent(in) :: local(:)
      real(wp), intent(in) :: centers(:, :)
      real(wp), intent(in) :: radii(:)
      real(wp), intent(in) :: reach
      integer, intent(in) :: nthreads
      type(seed_build_type), intent(out) :: geometry

      integer :: ic, r, n

      n = 0
      do ic = 1, size(local)
         do r = 1, local(ic)%batch%nrun
            if (local(ic)%seeded(r)) n = n + 1
         end do
      end do
      geometry%reach = reach
      allocate (geometry%anchors(3, n), geometry%radius_max(n))
      n = 0
      do ic = 1, size(local)
         do r = 1, local(ic)%batch%nrun
            if (.not. local(ic)%seeded(r)) cycle
            n = n + 1
            geometry%anchors(:, n) = local(ic)%batch%run_xyz(:, r)
            geometry%radius_max(n) = local(ic)%radius(r)
            if (local(ic)%grows(r)) geometry%radius_max(n) = max(local(ic)%radius(r), reach)
         end do
      end do
      geometry%centers = centers
      geometry%radii = radii
      geometry%nthreads = nthreads
   end subroutine search_geometry

   !> Order of the chunks in the search pass: estimated cost descending, ties by
   !> chunk index
   !>
   !> - cost: `extended_cost` per extended search plus one per window search
   !> - heaviest first keeps an expensive chunk from starting last and idling
   !>   the other threads; the rows do not depend on the order
   !>
   !> @param[in]  local Local passes of the chunks, in chunk order
   !> @param[out] order Permutation of `1:size(local)`
   subroutine search_order(local, order)
      type(local_pass_type), intent(in) :: local(:)
      integer, intent(out) :: order(:)

      integer, allocatable :: key(:)
      integer :: ic, r

      allocate (key(size(local)), source=0)
      do ic = 1, size(local)
         do r = 1, local(ic)%batch%nrun
            if (local(ic)%grows(r)) then
               key(ic) = key(ic) - extended_cost
            else if (local(ic)%seeded(r)) then
               key(ic) = key(ic) - 1
            end if
         end do
      end do
      call sort_index(key, order)
   end subroutine search_order

   !> Configure the worker from the cavity parameters
   !>
   !> - the LSF clone may borrow read-only data of `lsf_model`, which must outlive
   !>   the worker unchanged
   !> - the worker makes its own seeder from `param%seed` and a scratch of it
   !>
   !> @param[inout] self      Worker
   !> @param[in]    param     DROP parameters (validated)
   !> @param[in]    lsf_model Level set, already updated for the geometry; cloned
   !> @param[in]    centers   Atom centers (3, nat)
   !> @param[out]   error     Invalid level or seeder configuration
   subroutine worker_init(self, param, lsf_model, centers, error)
      class(projection_worker_type), intent(inout) :: self
      type(moist_cavity_drop_parameters_type), intent(in) :: param
      class(moist_cavity_drop_lsf_type), intent(in), target :: lsf_model
      real(wp), intent(in) :: centers(:, :)
      type(error_type), allocatable, intent(out) :: error

      if (.not. is_valid_proj_level(param%proj_level)) then
         call fatal_error(error, "DROP projection: level must be 1, 2, 3 or 4")
         return
      end if
      self%level = param%proj_level

      if (allocated(self%lsf)) deallocate (self%lsf)
      call lsf_model%clone_worker(self%lsf)
      call self%lsf%set_max_deriv(2)
      self%centers = centers
      self%reach_off = param%branch%reach_off
      self%full_reach = param%branch%full_reach

      self%solve = param%solve
      self%solve_seeded = param%solve
      self%solve_seeded%capture_duplicates = .true.
      self%polish = param%polish
      self%finalize = param%finalize
      self%phi_alpha = param%phi_alpha
      self%batch_size = param%batch_size

      call new_smoothstep_sigmoid(self%f_crit, param%w_0ls_off, param%w_0ls_on, &
                                  p_hi=param%w_0ls_p, a_hi=param%w_0ls_a, p_lo=param%w_0ls_p, &
                                  a_lo=param%w_0ls_a)
      call new_smoothstep_sigmoid(self%f_foc, param%w_0tra_off, param%w_0tra_on, &
                                  p_hi=param%w_0ls_p, a_hi=param%w_0ls_a, p_lo=param%w_0ls_p, &
                                  a_lo=param%w_0ls_a)
      call new_smoothstep_sigmoid(self%f_hor, param%hor_off, param%hor_on, &
                                  p_hi=param%w_0ls_p, a_hi=param%w_0ls_a, p_lo=param%w_0ls_p, &
                                  a_lo=param%w_0ls_a)
      ! Falling: open up to `reach_on`, closed from `reach_off`
      call new_smoothstep_sigmoid(self%f_reach, param%branch%reach_off, param%branch%reach_on, &
                                  p_hi=param%w_0ls_p, a_hi=param%w_0ls_a, p_lo=param%w_0ls_p, &
                                  a_lo=param%w_0ls_a)
      if (param%finalize%wleb_prune) then
         call new_smoothstep_sigmoid(self%f_wleb, param%finalize%wleb_prune_off, param%finalize%wleb_prune_on)
      else
         call new_smoothstep_sigmoid(self%f_wleb, 0.0_wp, 1.0_wp)
      end if
      call self%branch_weight%init(param%branch%dpdist_max, param%phi_alpha, param%branch%sep_cut)

      call self%seeder%init(param%seed, error)
      if (allocated(error)) return

      if (.not. allocated(self%passn)) allocate (self%passn)
      if (.not. allocated(self%passn_last)) allocate (self%passn_last)
      call self%passn%init(param%batch_size, error)
      if (allocated(error)) return
      call self%passn_last%init(param%batch_size, error)
      if (allocated(error)) return
      call self%pass2%init(param%batch_size, error)
      if (allocated(error)) return
      call self%rows%clear()
      self%n_seeded = 0
      self%n_seeds = 0
      self%n_skipped = 0
      self%n_escalated = 0
      self%n_captured = 0
      self%n_extended = 0
      self%n_window_skipped = 0
      self%pdist_search_max = 0.0_wp
      self%n_full = 0
   end subroutine worker_init

   !> Largest distance inside the branch window behind a fully open point
   !>
   !> @param[in] self  Worker
   !> @param[in] pdist Distance of the point (Bohr)
   !> @return    radius Window bound, at most the reach (Bohr)
   pure function worker_window_radius(self, pdist) result(radius)
      class(projection_worker_type), intent(in) :: self
      real(wp), intent(in) :: pdist
      real(wp) :: radius

      radius = min(self%reach_off, self%branch_weight%pdist_max(pdist))
   end function worker_window_radius

   !> Finalize record of one anchor
   !>
   !> @param[in]  self   Worker
   !> @param[in]  slot   Anchor slot, written into the rows
   !> @param[in]  anchor Anchor position
   !> @param[in]  owner  Owning atom
   !> @param[in]  wleb   Anchor quadrature weight
   !> @param[in]  f      Anchor switching value
   !> @param[out] rec    Per-anchor scalars
   pure subroutine anchor_record(self, slot, anchor, owner, wleb, f, rec)
      type(projection_worker_type), intent(in) :: self
      integer, intent(in) :: slot
      real(wp), intent(in) :: anchor(3)
      integer, intent(in) :: owner
      real(wp), intent(in) :: wleb, f
      type(finalize_anchor_type), intent(out) :: rec

      rec%index = slot
      rec%xyz = anchor
      rec%owner = owner
      rec%owner_center = self%centers(:, owner)
      rec%wleb = wleb
      rec%f = f
   end subroutine anchor_record

   !> Gate of an accepted column, the one the finalize will see ([[point_gate]])
   !>
   !> @param[in]  self  Worker
   !> @param[in]  batch Batch holding the column
   !> @param[in]  jcol  Accepted column
   !> @param[in]  rec   Finalize record of its anchor
   !> @param[out] gate  Gate
   !> @param[out] error Singular tangent Jacobian
   subroutine column_gate(self, batch, jcol, rec, gate, error)
      type(projection_worker_type), intent(in) :: self
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: jcol
      type(finalize_anchor_type), intent(in) :: rec
      real(wp), intent(out) :: gate
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: jac

      call point_gate(self%finalize, self%phi_alpha, batch, jcol, rec, self%f_crit, self%f_foc, self%f_hor, &
                      self%f_wleb, self%f_reach, jac, gate, error)
   end subroutine column_gate

   !* ================================================================================= *!
   !*                                     Local pass                                    *!
   !* ================================================================================= *!

   !> Local pass of one chunk of anchors
   !>
   !> - local Newton and polish from every anchor above the switching cut
   !> - per run, the accepted local column, the next step and the first search
   !>   radius
   !>
   !> @param[inout] self    Worker
   !> @param[in]    ctx     Run context, read only (polish traces)
   !> @param[in]    slots   Anchor slots, written into the rows (n)
   !> @param[in]    anchors Anchor positions (3, n)
   !> @param[in]    owners  Owning atoms (n)
   !> @param[in]    wlebs   Anchor quadrature weights (n)
   !> @param[in]    fs      Anchor switching values (n); below `wleb_cut` the anchor is skipped
   !> @param[inout] local   Local pass of the chunk, storage reused
   !> @param[out]   error   Build-fatal failure
   subroutine worker_solve_local(self, ctx, slots, anchors, owners, wlebs, fs, local, error)
      class(projection_worker_type), intent(inout) :: self
      type(moist_context_type), intent(in) :: ctx
      integer, intent(in) :: slots(:)
      real(wp), intent(in) :: anchors(:, :)
      integer, intent(in) :: owners(:)
      real(wp), intent(in) :: wlebs(:)
      real(wp), intent(in) :: fs(:)
      type(local_pass_type), intent(inout) :: local
      type(error_type), allocatable, intent(out) :: error

      type(finalize_anchor_type) :: rec
      integer :: k, r, j, jacc, icol
      real(wp) :: gate, window

      ! A chunk's results must not depend on the chunks this thread solved before
      call self%lsf%reset_history()
      if (local%batch%capacity == 0) then
         call local%batch%init(self%batch_size, error)
         if (allocated(error)) return
      else
         call local%batch%clear()
      end if
      do k = 1, size(slots)
         if (fs(k) < self%finalize%wleb_cut) cycle
         call local%batch%begin_run(k, slots(k), anchors(:, k), owners(k))
         call local%batch%append(anchors(:, k), 0.0_wp, icol)
      end do
      if (local%batch%n == 0) return
      call newton_solve_batch(local%batch, self%lsf, self%solve, self%polish, self%phi_alpha, self%nscratch, error)
      if (allocated(error)) return
      call polish_batch(ctx, local%batch, self%lsf, self%polish, self%solve, self%phi_alpha, &
                        self%nscratch, error)
      if (allocated(error)) return

      if (allocated(local%next)) then
         if (size(local%next) < local%batch%nrun) deallocate (local%next, local%jloc, local%radius)
      end if
      if (.not. allocated(local%next)) then
         allocate (local%next(local%batch%nrun), local%jloc(local%batch%nrun), local%radius(local%batch%nrun))
      end if

      ! Only a fully open local point bounds the search by its window
      do r = 1, local%batch%nrun
         k = local%batch%run_key(r)
         j = local%batch%run_start(r)
         jacc = 0
         if (local%batch%status(j) == point_accepted) jacc = j
         local%jloc(r) = jacc
         if (jacc == 0) then
            local%next(r) = next_extended
            local%radius(r) = first_ball(r, 0)
         else if (self%level == 1) then
            local%next(r) = next_finalize
            local%radius(r) = 0.0_wp
         else
            call anchor_record(self, slots(k), anchors(:, k), owners(k), wlebs(k), fs(k), rec)
            call column_gate(self, local%batch, jacc, rec, gate, error)
            if (allocated(error)) return
            if (gate == 1.0_wp .and. .not. self%full_reach) then
               window = self%window_radius(local%batch%pdist(jacc))
               if (window <= 0.0_wp) then
                  ! Nothing behind a zero window carries weight
                  local%next(r) = next_finalize
                  local%radius(r) = 0.0_wp
               else if (self%level == 2 .and. local%batch%pdist(jacc) <= level2_window_pdist) then
                  local%next(r) = next_finalize
                  local%radius(r) = 0.0_wp
                  self%n_window_skipped = self%n_window_skipped + 1
               else
                  local%next(r) = next_window
                  local%radius(r) = window
               end if
            else
               local%next(r) = next_extended
               local%radius(r) = first_ball(r, jacc)
            end if
         end if
      end do

   contains

      !> First ball of an extended search: the window ball of the local point,
      !> else of the anchor's Newton length `|S|/|grad S|`, a guess the search
      !> outgrows when wrong
      !>
      !> @param[in] r1   Run index
      !> @param[in] jloc Accepted local column, 0 for none
      !> @return    ball Radius (Bohr)
      function first_ball(r1, jloc) result(ball)
         integer, intent(in) :: r1, jloc
         real(wp) :: ball

         real(wp) :: s_anchor, g_anchor(3)

         if (jloc > 0) then
            ball = self%window_radius(local%batch%pdist(jloc))
         else
            call local%batch%anchor_jet(r1, s_anchor, g_anchor)
            ball = self%reach_off
            if (norm2(g_anchor) > 0.0_wp) ball = self%window_radius(min(self%reach_off, abs(s_anchor)/norm2(g_anchor)))
         end if
         if (self%full_reach) ball = self%reach_off
         ball = max(ball, min_ball_fraction*self%reach_off)
      end function first_ball

   end subroutine worker_solve_local

   !* ================================================================================= *!
   !*                                    Search pass                                    *!
   !* ================================================================================= *!

   !> Search pass of one chunk of anchors, from its local pass
   !>
   !> - `seeder` may be the worker's own `self%seeder`: the pass must never
   !>   define that component, as the two dummies then name one object
   !>
   !> @param[inout] self     Worker; its `scratch` takes the seeds
   !> @param[in]    ctx      Run context, read only (polish traces)
   !> @param[in]    local    Local pass of the chunk
   !> @param[in]    seeder   Built seeder, read only
   !> @param[in]    slots    Anchor slots, written into the rows (n)
   !> @param[in]    anchors  Anchor positions (3, n)
   !> @param[in]    owners   Owning atoms (n)
   !> @param[in]    wlebs    Anchor quadrature weights (n)
   !> @param[in]    fs       Anchor switching values (n)
   !> @param[out]   nrows    Rows appended per anchor (n)
   !> @param[out]   lost     Anchor kept no accepted point (n)
   !> @param[out]   error    Build-fatal failure
   !> @param[out]   unseeded Lost anchor whose searches up to the reach gave no seed (n)
   subroutine worker_search(self, ctx, local, seeder, slots, anchors, owners, wlebs, fs, nrows, lost, error, &
                            unseeded)
      class(projection_worker_type), intent(inout) :: self
      type(moist_context_type), intent(in) :: ctx
      type(local_pass_type), intent(in) :: local
      type(projection_seeder_type), intent(in) :: seeder
      integer, intent(in) :: slots(:)
      real(wp), intent(in) :: anchors(:, :)
      integer, intent(in) :: owners(:)
      real(wp), intent(in) :: wlebs(:)
      real(wp), intent(in) :: fs(:)
      integer, intent(out) :: nrows(:)
      logical, intent(out) :: lost(:)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(out), optional :: unseeded(:)

      integer :: k, r, j, icol, nnear, nnear_seeds

      ! A chunk's results must not depend on the chunks this thread searched before
      call self%lsf%reset_history()
      nrows = 0
      lost = .false.
      if (present(unseeded)) unseeded = .false.
      if (local%batch%n == 0) return

      call self%pass2%clear()
      call self%passn%clear()
      nnear = 0
      nnear_seeds = 0
      do r = 1, local%batch%nrun
         k = local%batch%run_key(r)
         j = local%batch%run_start(r)
         select case (local%next(r))
         case (next_finalize)
            call finish_anchor(local%batch, j, j, k)
         case (next_window)
            call seed_anchor(local%jloc(r), k, local%radius(r))
         case default
            if (local%jloc(r) == 0) self%n_escalated = self%n_escalated + 1
            call queue_search(r, k, local%jloc(r), local%radius(r))
         end select
         if (allocated(error)) return
      end do
      call solve_queue()
      if (allocated(error)) return
      call flush_pass2()

   contains

      !> Request of anchor `k1`: the whole ball of radius `radius`
      !>
      !> - the radius is positive: a window ball of a fully open point (a zero
      !>   window finalizes in the local pass) or an extended-search ball, at
      !>   least `min_ball_fraction` of the reach
      !> - the seed filters of the level are the seeder's settings
      !>
      !> @param[in]  k1           Anchor index within the chunk
      !> @param[in]  radius       Search radius (Bohr), positive
      !> @param[in]  has_accepted Whether the anchor already has an accepted local point
      !> @param[out] request      Request
      subroutine anchor_request(k1, radius, has_accepted, request)
         integer, intent(in) :: k1
         real(wp), intent(in) :: radius
         logical, intent(in) :: has_accepted
         type(seed_request_type), intent(out) :: request

         request%anchor = anchors(:, k1)
         request%owner = owners(k1)
         request%owner_center = self%centers(:, owners(k1))
         request%pdist_start = radius
         request%has_accepted = has_accepted
      end subroutine anchor_request

      !> Seed the window ball of a fully open local point and open its pass-2 run
      !>
      !> - flushes the batch first when full
      !>
      !> @param[in] jacc   Local column of the accepted local point
      !> @param[in] k1     Anchor index within the chunk
      !> @param[in] radius Window radius behind the local point (Bohr)
      subroutine seed_anchor(jacc, k1, radius)
         !> Local column of the accepted local point
         integer, intent(in) :: jacc
         !> Anchor index within the chunk
         integer, intent(in) :: k1
         !> Window radius behind the local point (Bohr)
         real(wp), intent(in) :: radius
         type(seed_request_type) :: request
         integer :: n_new, i

         self%pdist_search_max = max(self%pdist_search_max, radius)
         call anchor_request(k1, radius, .true., request)
         call seeder%seed(self%scratch, self%lsf, request, error)
         if (allocated(error)) return
         self%n_seeded = self%n_seeded + 1
         self%n_seeds = self%n_seeds + self%scratch%n_points

         n_new = self%scratch%n_points + 1
         if (self%pass2%would_overflow(n_new)) then
            call flush_pass2()
            if (allocated(error)) return
         end if
         call self%pass2%begin_run(k1, slots(k1), anchors(:, k1), owners(k1))
         call self%pass2%append_copy(local%batch, jacc, icol)
         ! Seeds start with a zero multiplier
         do i = 1, self%scratch%n_points
            call self%pass2%append(self%scratch%points(:, i), 0.0_wp, icol)
         end do
      end subroutine seed_anchor

      !> Queue the anchor of local run `r1` for the extended search
      !>
      !> @param[in] r1   Local run index
      !> @param[in] k1   Anchor index within the chunk
      !> @param[in] jacc Local column of the accepted local point, 0 for none
      !> @param[in] ball First ball (Bohr)
      subroutine queue_search(r1, k1, jacc, ball)
         integer, intent(in) :: r1, k1, jacc
         real(wp), intent(in) :: ball

         if (ball >= self%reach_off) self%n_full = self%n_full + 1

         nnear = nnear + 1
         call reserve_queue(nnear, nnear_seeds)
         ! Nothing solved yet
         self%queue(nnear) = queue_entry_type(run=r1, jloc=jacc, first=nnear_seeds + 1, next=nnear_seeds + 1, &
                                              take=self%nearest_count, radius=self%reach_off, ball=ball)
         self%n_seeded = self%n_seeded + 1
         self%n_extended = self%n_extended + 1
         call seed_ball(nnear, k1, -1.0_wp)
      end subroutine queue_search

      !> Seed the ball (or outer shell) of a queued anchor, growing it while empty
      !>
      !> - the seeds replace the anchor's unsolved ones
      !> - when the seeder supports shell queries (`scratch%shell_filtered`), seeds
      !>   strictly outside `prev_ball` are returned and `mark_solved` is skipped
      !> - otherwise a seed equal to one solved from the anchor's previous list is
      !>   marked done: a grown ball proposes the seeds of the smaller one again
      !>
      !> @param[in] e         Queue index of the anchor
      !> @param[in] k1        Anchor index within the chunk
      !> @param[in] prev_ball Inner ball already exhausted (Bohr), negative for full ball
      subroutine seed_ball(e, k1, prev_ball)
         integer, intent(in) :: e, k1
         real(wp), intent(in) :: prev_ball
         type(seed_request_type) :: request
         integer :: n, solved_first, solved_last
         real(wp) :: inner

         ! Solved so far: the head of the current list
         solved_first = self%queue(e)%first
         solved_last = self%queue(e)%next - 1
         inner = prev_ball
         do
            call anchor_request(k1, self%queue(e)%ball, self%queue(e)%jloc > 0, request)
            request%nearest_first = .true.
            request%pdist_inner = inner
            call seeder%seed(self%scratch, self%lsf, request, error)
            if (allocated(error)) return
            n = self%scratch%n_points
            call reserve_queue(nnear, nnear_seeds + n)
            ! `points` stays unallocated until the first seed
            if (n > 0) self%near_seeds(:, nnear_seeds + 1:nnear_seeds + n) = self%scratch%points(:, 1:n)
            self%near_done(nnear_seeds + 1:nnear_seeds + n) = .false.
            if (.not. self%scratch%shell_filtered) then
               call mark_solved(solved_first, solved_last, nnear_seeds + 1, nnear_seeds + n)
            end if
            self%queue(e)%reach = self%scratch%reach
            self%queue(e)%first = nnear_seeds + 1
            self%queue(e)%next = nnear_seeds + 1
            nnear_seeds = nnear_seeds + n
            self%queue(e)%last = nnear_seeds
            if (n > 0 .or. self%queue(e)%ball >= self%queue(e)%radius) exit
            inner = self%queue(e)%ball
            self%queue(e)%ball = min(self%queue(e)%radius, 2.0_wp*self%queue(e)%ball)
         end do
      end subroutine seed_ball

      !> Mark the seeds `first:last` equal to a solved seed `solved_first:solved_last` as done
      !>
      !> - exact equality: only a start solved before is skipped
      !> - solved seeds sorted by their x coordinate, an exact key; points of an
      !>   equal key compared whole
      !>
      !> @param[in] solved_first First solved seed
      !> @param[in] solved_last  Last solved seed
      !> @param[in] first        First new seed
      !> @param[in] last         Last new seed
      subroutine mark_solved(solved_first, solved_last, first, last)
         integer, intent(in) :: solved_first, solved_last, first, last
         real(wp), allocatable :: key(:)
         integer, allocatable :: order(:)
         integer :: nsolved, s, lo, hi, mid, i

         nsolved = solved_last - solved_first + 1
         if (nsolved <= 0 .or. last < first) return
         key = self%near_seeds(1, solved_first:solved_last)
         allocate (order(nsolved))
         call sort_index(key, order)
         do s = first, last
            ! First solved key not below the seed's
            lo = 1
            hi = nsolved + 1
            do while (lo < hi)
               mid = (lo + hi)/2
               if (key(mid) < self%near_seeds(1, s)) then
                  lo = mid + 1
               else
                  hi = mid
               end if
            end do
            do i = lo, nsolved
               if (key(i) > self%near_seeds(1, s)) exit
               if (all(self%near_seeds(:, solved_first + order(i) - 1) == self%near_seeds(:, s))) then
                  self%near_done(s) = .true.
                  exit
               end if
            end do
         end do
      end subroutine mark_solved

      !> Grow the search queue to `n` anchors and `m` seeds, keeping its contents
      !>
      !> @param[in] n Anchors to hold
      !> @param[in] m Seeds to hold
      subroutine reserve_queue(n, m)
         integer, intent(in) :: n, m
         type(queue_entry_type), allocatable :: tmp(:)
         integer :: cap

         if (.not. allocated(self%queue)) then
            allocate (self%queue(64), self%near_seeds(3, 1024), self%near_done(1024))
         end if
         if (n > size(self%queue)) then
            allocate (tmp(max(n, 2*size(self%queue))))
            tmp(1:size(self%queue)) = self%queue
            call move_alloc(tmp, self%queue)
         end if
         if (m > size(self%near_seeds, 2)) then
            cap = max(m, 2*size(self%near_seeds, 2))
            call grow_array(self%near_seeds, 3, cap)
            call grow_array(self%near_done, cap)
         end if
      end subroutine reserve_queue

      !> Extended search of every queued anchor, in rounds batched across anchors
      !>
      !> - each pending anchor takes `take` seeds (starting at `nearest_count`,
      !>   doubling per active round, reset on ball growth at level 1) that lie
      !>   within the seeder's reach of its radius, skipping solved ones
      !> - after the round the radius of an anchor shrinks to the window behind
      !>   its nearest qualifying accepted point (`gate == 1` at levels 2-4,
      !>   `gate > 0` at level 1)
      !> - an anchor with no unsolved seed left inside gets a larger ball while
      !>   its ball is short of its radius, and else is finalized from its
      !>   latest run, which holds every point accepted for it
      subroutine solve_queue()
         type(finalize_anchor_type) :: rec
         integer :: npend, nkeep, ip, e, k1, rn, s, s_last, jn, j_first, j_last, col0, run0, nskip
         integer, allocatable :: pend(:), run_new_first(:)
         real(wp) :: limit, gate_j, prev_ball
         logical :: more

         if (nnear == 0) return
         ! A round without seeds would leave every waiting anchor waiting forever
         if (self%nearest_count < 1) then
            call fatal_error(error, "DROP projection: nearest_count must be positive")
            return
         end if
         allocate (pend(nnear), run_new_first(nnear))
         do e = 1, nnear
            pend(e) = e
         end do
         call self%ext_saddles%reset(nnear)
         npend = nnear
         do while (npend > 0)
            call carry_latest_runs(pend(1:npend))
            col0 = self%passn%n + 1
            run0 = self%passn%nrun + 1
            do ip = 1, npend
               e = pend(ip)
               k1 = local%batch%run_key(self%queue(e)%run)
               limit = self%queue(e)%radius + self%queue(e)%reach
               s_last = self%queue(e)%next - 1
               do s = self%queue(e)%next, min(self%queue(e)%last, self%queue(e)%next + self%queue(e)%take - 1)
                  if (norm2(self%near_seeds(:, s) - anchors(:, k1)) > limit) exit
                  s_last = s
               end do
               if (s_last < self%queue(e)%next) cycle
               nskip = count(self%near_done(self%queue(e)%next:s_last))
               self%n_skipped = self%n_skipped + nskip
               ! Every seed of the round solved before: no run, the latest one holds the accepted points
               if (nskip == s_last - self%queue(e)%next + 1) then
                  self%queue(e)%next = s_last + 1
                  if (self%queue(e)%take <= max_near_take) self%queue(e)%take = 2*self%queue(e)%take
                  cycle
               end if

               ! Points accepted so far lead the run: kept, and captors of the new seeds
               j_first = 1
               j_last = 0
               if (self%queue(e)%tail > 0) then
                  j_first = self%passn%run_start(self%queue(e)%tail)
                  j_last = j_first + self%passn%run_len(self%queue(e)%tail) - 1
               end if
               call self%passn%begin_run(e, slots(k1), anchors(:, k1), owners(k1))
               if (self%queue(e)%tail > 0) then
                  do jn = j_first, j_last
                     if (self%passn%status(jn) == point_accepted) call self%passn%append_own(jn, icol)
                  end do
               else if (self%queue(e)%jloc > 0) then
                  call self%passn%append_copy(local%batch, self%queue(e)%jloc, icol)
               end if
               self%queue(e)%tail = self%passn%nrun
               run_new_first(self%passn%nrun - run0 + 1) = self%passn%n + 1
               do s = self%queue(e)%next, s_last
                  if (.not. self%near_done(s)) call self%passn%append(self%near_seeds(:, s), 0.0_wp, icol)
               end do
               self%n_seeds = self%n_seeds + s_last - self%queue(e)%next + 1 - nskip
               self%queue(e)%next = s_last + 1
               if (self%queue(e)%take <= max_near_take) self%queue(e)%take = 2*self%queue(e)%take
            end do

            if (self%passn%n >= col0) then
               call newton_solve_batch(self%passn, self%lsf, self%solve_seeded, self%polish, self%phi_alpha, &
                                       self%nscratch, error, first=col0)
               if (allocated(error)) return
               call polish_batch(ctx, self%passn, self%lsf, self%polish, self%solve_seeded, self%phi_alpha, &
                                 self%nscratch, error, first=col0, saddle_cache=self%ext_saddles)
               if (allocated(error)) return
               self%n_captured = self%n_captured + count(self%passn%status(col0:self%passn%n) == point_captured)
            end if

            ! Radius behind the nearest qualifying newly accepted point:
            ! copied columns from prior rounds and the local point were already checked
            if (.not. self%full_reach) then
               do rn = run0, self%passn%nrun
                  e = self%passn%run_key(rn)
                  k1 = local%batch%run_key(self%queue(e)%run)
                  call anchor_record(self, slots(k1), anchors(:, k1), owners(k1), wlebs(k1), fs(k1), rec)
                  do jn = run_new_first(rn - run0 + 1), self%passn%run_start(rn) + self%passn%run_len(rn) - 1
                     if (self%passn%status(jn) /= point_accepted) cycle
                     call column_gate(self, self%passn, jn, rec, gate_j, error)
                     if (allocated(error)) return
                     if (gate_j == 1.0_wp .or. (self%level == 1 .and. gate_j > 0.0_wp)) then
                        self%queue(e)%radius = min(self%queue(e)%radius, self%window_radius(self%passn%pdist(jn)))
                     end if
                  end do
               end do
            end if

            nkeep = 0
            do ip = 1, npend
               e = pend(ip)
               k1 = local%batch%run_key(self%queue(e)%run)
               more = self%queue(e)%next <= self%queue(e)%last
               if (more) then
                  more = norm2(self%near_seeds(:, self%queue(e)%next) - anchors(:, k1)) &
                         <= self%queue(e)%radius + self%queue(e)%reach
               end if
               if (more) then
                  nkeep = nkeep + 1
                  pend(nkeep) = e
               else if (self%queue(e)%ball < self%queue(e)%radius) then
                  ! To the radius a qualifying point set, else twice as far
                  prev_ball = self%queue(e)%ball
                  if (self%queue(e)%next <= self%queue(e)%last) prev_ball = -1.0_wp
                  if (self%queue(e)%radius < self%reach_off) then
                     self%queue(e)%ball = self%queue(e)%radius
                  else
                     self%queue(e)%ball = min(self%queue(e)%radius, 2.0_wp*self%queue(e)%ball)
                  end if
                  if (self%level == 1) self%queue(e)%take = self%nearest_count
                  call seed_ball(e, k1, prev_ball)
                  nkeep = nkeep + 1
                  pend(nkeep) = e
               else
                  call finish_queued(e, k1)
               end if
               if (allocated(error)) return
            end do
            npend = nkeep
         end do
      end subroutine solve_queue

      !> Start a round's batch with the latest run of every pending anchor,
      !> reduced to its accepted points, from the last round's batch
      !>
      !> - superseded runs, finalized anchors and columns not accepted are dropped,
      !>   so the batch holds the live runs only; the new runs follow them
      !> - the finalize and the next round read only the accepted points of a run
      !>
      !> @param[in] pend Pending anchors, queue indices
      subroutine carry_latest_runs(pend)
         integer, intent(in) :: pend(:)
         type(projection_batch_type), allocatable :: last
         integer :: ip, e, rn, jn

         call move_alloc(self%passn, last)
         call move_alloc(self%passn_last, self%passn)
         call self%passn%clear()
         do ip = 1, size(pend)
            e = pend(ip)
            rn = self%queue(e)%tail
            if (rn == 0) cycle
            call self%passn%begin_run(e, last%run_id(rn), last%run_xyz(:, rn), last%run_owner(rn))
            do jn = last%run_start(rn), last%run_start(rn) + last%run_len(rn) - 1
               if (last%status(jn) == point_accepted) call self%passn%append_copy(last, jn, icol)
            end do
            self%queue(e)%tail = self%passn%nrun
         end do
         call move_alloc(last, self%passn_last)
      end subroutine carry_latest_runs

      !> Finalize a queued anchor from its latest run, or from its local run without one
      !>
      !> @param[in] e  Queue index of the anchor
      !> @param[in] k1 Anchor index within the chunk
      subroutine finish_queued(e, k1)
         integer, intent(in) :: e, k1
         integer :: rn, j1

         self%pdist_search_max = max(self%pdist_search_max, self%queue(e)%ball)
         if (self%queue(e)%tail > 0) then
            rn = self%queue(e)%tail
            call finish_anchor(self%passn, self%passn%run_start(rn), &
                               self%passn%run_start(rn) + self%passn%run_len(rn) - 1, k1)
         else
            j1 = local%batch%run_start(self%queue(e)%run)
            call finish_anchor(local%batch, j1, j1, k1)
            ! No seed was ever solved: no sign change within the reach
            if (present(unseeded)) unseeded(k1) = lost(k1)
         end if
      end subroutine finish_queued

      !> Solve, polish and finalize every run of the pass-2 batch
      subroutine flush_pass2()
         integer :: r2, k2, first, last

         if (self%pass2%nrun == 0) return
         if (self%pass2%n > 0) then
            call newton_solve_batch(self%pass2, self%lsf, self%solve_seeded, self%polish, self%phi_alpha, &
                                    self%nscratch, error)
            if (allocated(error)) return
            call polish_batch(ctx, self%pass2, self%lsf, self%polish, self%solve_seeded, self%phi_alpha, &
                              self%nscratch, error)
            if (allocated(error)) return
            self%n_captured = self%n_captured + count(self%pass2%status(1:self%pass2%n) == point_captured)
         end if
         do r2 = 1, self%pass2%nrun
            k2 = self%pass2%run_key(r2)
            first = self%pass2%run_start(r2)
            last = first + self%pass2%run_len(r2) - 1
            call finish_anchor(self%pass2, first, last, k2)
            if (allocated(error)) return
         end do
         call self%pass2%clear()
      end subroutine flush_pass2

      !> Finalize one anchor's column range into rows
      !>
      !> @param[in] batch Batch holding the anchor's columns
      !> @param[in] first First column of the anchor's run
      !> @param[in] last  Last column of the anchor's run
      !> @param[in] k2    Anchor index within the chunk
      subroutine finish_anchor(batch, first, last, k2)
         type(projection_batch_type), intent(in) :: batch
         integer, intent(in) :: first, last
         integer, intent(in) :: k2
         type(finalize_anchor_type) :: rec
         integer :: npoint

         call anchor_record(self, slots(k2), anchors(:, k2), owners(k2), wlebs(k2), fs(k2), rec)
         call finalize_anchor(self%finalize, self%phi_alpha, batch, first, last, rec, self%f_crit, &
                              self%f_foc, self%f_hor, self%f_wleb, self%f_reach, self%branch_weight, self%rows, &
                              npoint, nrows(k2), self%fscratch, error)
         lost(k2) = npoint == 0
      end subroutine finish_anchor

   end subroutine worker_search

   !> Run both passes on one chunk of anchors
   !>
   !> - between the passes the worker's own seeder is built from the chunk's
   !>   search geometry ([[search_geometry]]), unless `seeder` names a built one
   !> - the scratch is made anew from the seeder in use
   !>
   !> @param[inout] self     Worker
   !> @param[in]    ctx      Run context, read only (polish traces, build threads)
   !> @param[in]    slots    Anchor slots, written into the rows (n)
   !> @param[in]    anchors  Anchor positions (3, n)
   !> @param[in]    owners   Owning atoms (n)
   !> @param[in]    wlebs    Anchor quadrature weights (n)
   !> @param[in]    fs       Anchor switching values (n); below `wleb_cut` the anchor is skipped
   !> @param[out]   nrows    Rows appended per anchor (n)
   !> @param[out]   lost     Anchor kept no accepted point (n)
   !> @param[out]   error    Build-fatal failure
   !> @param[out]   unseeded Lost anchor whose searches up to the reach gave no seed (n)
   !> @param[in]    seeder   Built seeder to search with instead of the worker's own
   subroutine worker_process_chunk(self, ctx, slots, anchors, owners, wlebs, fs, nrows, lost, error, unseeded, &
                                   seeder)
      class(projection_worker_type), intent(inout) :: self
      type(moist_context_type), intent(in) :: ctx
      integer, intent(in) :: slots(:)
      real(wp), intent(in) :: anchors(:, :)
      integer, intent(in) :: owners(:)
      real(wp), intent(in) :: wlebs(:)
      real(wp), intent(in) :: fs(:)
      integer, intent(out) :: nrows(:)
      logical, intent(out) :: lost(:)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(out), optional :: unseeded(:)
      type(projection_seeder_type), intent(in), optional :: seeder

      type(local_pass_type) :: local(1)
      type(seed_build_type) :: geometry

      nrows = 0
      lost = .false.
      if (present(unseeded)) unseeded = .false.
      call self%solve_local(ctx, slots, anchors, owners, wlebs, fs, local(1), error)
      if (allocated(error)) return
      if (present(seeder)) then
         call self%search(ctx, local(1), seeder, slots, anchors, owners, wlebs, fs, nrows, lost, error, unseeded)
         return
      end if
      call search_geometry(local, self%centers, self%lsf%radii, self%reach_off, ctx%get_num_threads(), geometry)
      if (size(geometry%radius_max) > 0) then
         call self%seeder%build(self%lsf, geometry, error)
         if (allocated(error)) return
      end if
      call self%search(ctx, local(1), self%seeder, slots, anchors, owners, wlebs, fs, nrows, lost, error, unseeded)
   end subroutine worker_process_chunk

   !> Project one anchor and return its rows (the single-anchor entry for tests)
   !>
   !> - the worker's row buffer is cleared first, so the outputs are exactly
   !>   this anchor's surviving branches in canonical order
   !> - the worker's own seeder, built for this anchor
   !> - without rows the row outputs are allocated empty; on error none is allocated
   !>
   !> @param[inout] self      Worker
   !> @param[in]    ctx       Run context, read only (polish traces)
   !> @param[in]    anchor    Anchor position
   !> @param[in]    owner     Owning atom
   !> @param[out]   n_points  Rows produced (0: anchor lost or below the cut)
   !> @param[out]   error     Build-fatal failure
   !> @param[out]   points    Projected points (3, n_points)
   !> @param[out]   normals   Unit normals (3, n_points)
   !> @param[out]   lambda    Multipliers
   !> @param[out]   pdist     Anchor distances
   !> @param[out]   phi       Objective values
   !> @param[out]   weights   Branch weights
   !> @param[out]   n_seeds   Seeds solved for this anchor
   !> @param[out]   lost      No accepted point at all
   !> @param[in]    wleb      Anchor quadrature weight (default 1)
   !> @param[in]    f         Anchor switching value (default 1)
   !> @param[out]   unseeded  Lost, and no search gave a seed to solve
   subroutine worker_project_anchor(self, ctx, anchor, owner, n_points, error, points, normals, lambda, pdist, &
                                    phi, weights, n_seeds, lost, wleb, f, unseeded)
      class(projection_worker_type), intent(inout) :: self
      type(moist_context_type), intent(in) :: ctx
      real(wp), intent(in) :: anchor(3)
      integer, intent(in) :: owner
      integer, intent(out) :: n_points
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable, intent(out), optional :: points(:, :), normals(:, :)
      real(wp), allocatable, intent(out), optional :: lambda(:), pdist(:), phi(:), weights(:)
      integer, intent(out), optional :: n_seeds
      logical, intent(out), optional :: lost
      real(wp), intent(in), optional :: wleb
      real(wp), intent(in), optional :: f
      logical, intent(out), optional :: unseeded

      real(wp) :: wleb_k(1), f_k(1)
      integer :: nrows(1), seeds_before
      logical :: lost_k(1), unseeded_k(1)

      wleb_k = 1.0_wp
      f_k = 1.0_wp
      if (present(wleb)) wleb_k = wleb
      if (present(f)) f_k = f

      call self%rows%clear()
      seeds_before = self%n_seeds
      call self%process_chunk(ctx, [1], reshape(anchor, [3, 1]), [owner], wleb_k, f_k, nrows, lost_k, error, &
                              unseeded=unseeded_k)
      n_points = 0
      if (allocated(error)) return
      n_points = self%rows%n
      if (present(n_seeds)) n_seeds = self%n_seeds - seeds_before
      if (present(lost)) lost = lost_k(1)
      if (present(unseeded)) unseeded = unseeded_k(1)
      ! A fresh buffer without rows is unallocated: the empty sections need storage
      call self%rows%reserve(1)
      if (present(points)) points = self%rows%xyz(:, 1:n_points)
      if (present(normals)) normals = self%rows%normal(:, 1:n_points)
      if (present(lambda)) lambda = self%rows%lambda(1:n_points)
      if (present(pdist)) pdist = self%rows%pdist(1:n_points)
      if (present(phi)) phi = self%rows%phi(1:n_points)
      if (present(weights)) weights = self%rows%wbranch(1:n_points)
   end subroutine worker_project_anchor

end module moist_cavity_drop_projection_worker
