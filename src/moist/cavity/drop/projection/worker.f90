!> Per-thread worker of the batched DROP projection
!>
!> One worker runs the whole per-anchor pipeline on a chunk of anchors:
!>
!>   1. local Newton from every anchor (pass-1 batch), optimality polish; an
!>      anchor so deep inside another sphere that its path gates are closed
!>      for every point is left out
!>   2. gate of the accepted local point ([[point_gate]]), which sets how far
!>      the anchor is searched (levels 3, 4; level 1 searches no anchor with a
!>      local point)
!>   3. fully open local point: the seeder covers the ball up to the branch
!>      window behind it, at most the reach; its seeds, plus the local point,
!>      form one run in the pass-2 batch, flushed when full
!>   4. any other anchor: extended search (below)
!>   5. batched Newton and polish on the seeds; a seed whose Newton prediction
!>      lands on an accepted point of its run stops as a duplicate
!>   6. per-anchor finalize (dedupe, canonical order, CP Jacobian, gates,
!>      branch weights, cut) appending rows to the worker's row buffer
!>
!> Search radius of an anchor:
!>
!>   - a point beyond the reach `branch%reach_off` has a closed gate, and one
!>     more than the branch window behind a fully open point (gate exactly one)
!>     has no weight and gives none, so the radius is
!>     `min(reach, window(pdist_k))`, `k` the nearest fully open point found,
!>     and the reach when there is none
!>   - a point with a partly open or closed gate never shortens it: what lies
!>     behind keeps weight
!>
!> Extended search of an anchor whose local point is not fully open or missing:
!>
!>   - the seeder (level 1: the fallback) covers a ball around the anchor,
!>     seeds nearest first: the window ball of the local point, or of the
!>     anchor's Newton length without one
!>   - rounds batched over all such anchors of the chunk: the next seeds of
!>     every anchor, twice as many per round, are solved and the radius is
!>     taken again from the fully open points accepted so far
!>   - an anchor whose seeds within the seeder's `reach` of its radius are
!>     solved is done when its ball covers the radius; else the ball grows, to
!>     the radius when a fully open point set one, to twice its size (at most
!>     the reach) when none did, and is seeded again
!>   - every round's run of an anchor starts with copies of its points accepted
!>     so far, so the last run holds all of them and is the one finalized
!>
!> Levels (TODO: will be changed soon):
!>
!> | level | searched anchors | seeder        | failed local solve |
!> |-------|------------------|---------------|--------------------|
!> | 1     | none             | none          | coarse octree      |
!> | 3     | all              | coarse octree | coarse octree      |
!> | 4     | all              | fine octree    | fine octree         |
!>
!> - any other level is refused
!> - an anchor without an accepted local point goes to the `fallback` seeder
!>   (level 1: level-3 octree settings), else to the level's seeder; with
!>   seeds but nothing accepted its seeds are solved once more with the full
!>   system from the start (the cheap step can carry a seed off a thin neck);
!>   with nothing accepted at the end it is lost, reported through `lost`, and
!>   appends no row; a lost anchor whose searches gave no seed to solve has
!>   no surface within the reach and is also reported through `unseeded`
!> - a surviving leaf without sign-change evidence is dropped by the coarse
!>   octree and refined once by the fine one
!> - `nearest_count` is a starting guess, not a correctness parameter
!> - `full_reach` (debug) sends every searched anchor through the extended
!>   search with ball and radius held at the reach; the rows are those of the
!>   adaptive radius as far as the seeder finds the same points
!> - LSF failures, an exhausted octree budget, an anchor on its owner center
!>   and a singular tangent block are build-fatal and come back as `error`
module moist_cavity_drop_projection_worker
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type, is_valid_proj_level
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_projection_batch, only: projection_batch_type, point_accepted, point_captured
   use moist_cavity_drop_projection_parameters, only: drop_seed_parameters_type, drop_solve_parameters_type, &
                                                     drop_polish_parameters_type, drop_finalize_parameters_type
   use moist_cavity_drop_projection_newton, only: newton_scratch_type, newton_solve_batch
   use moist_cavity_drop_projection_polish, only: polish_batch
   use moist_cavity_drop_projection_finalize, only: finalize_anchor_type, &
                                                   finalize_scratch_type, projection_rows_type, &
                                                   finalize_anchor, point_gate
   use moist_cavity_drop_projection_seed, only: projection_seeder_type, seed_request_type
   use moist_cavity_drop_projection_seed_octree, only: octree_seeder_type
   implicit none(type, external)
   private

   public :: projection_worker_type

   !> Seeds per anchor in the first round of the extended search
   integer, parameter :: default_nearest_count = 4
   !> Smallest first ball of an extended search, in units of the reach
   real(wp), parameter :: min_ball_fraction = 0.02_wp

   !> Per-thread projection worker
   type :: projection_worker_type
      !> Projection level (1, 3 or 4)
      integer :: level = 3
      !> Thread-local level set clone
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      !> Atom centers (3, nat)
      real(wp), allocatable :: centers(:, :)
      !> Stage parameters
      type(drop_solve_parameters_type) :: solve
      type(drop_polish_parameters_type) :: polish
      type(drop_finalize_parameters_type) :: finalize
      !> Solve parameters of the seeded pass: `solve` with duplicate capture
      type(drop_solve_parameters_type) :: solve_seeded
      !> Solve parameters of the restart: `solve_seeded` without cheap steps
      type(drop_solve_parameters_type) :: solve_restart
      !> Anchor objective weight, set by `init`
      real(wp) :: phi_alpha
      !> Distance from which no point carries weight (Bohr): the largest search radius
      real(wp) :: reach_off = 0.0_wp
      !> Search every anchor to the reach, whatever its points (debug)
      logical :: full_reach = .false.
      !> Seeder of the level, searching every anchor; unallocated at level 1
      class(projection_seeder_type), allocatable :: seeder
      !> Seeder of anchors without an accepted local point; unallocated: `seeder`
      class(projection_seeder_type), allocatable :: fallback
      !> Switching functions and branch weight model
      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc, f_hor, f_wleb
      !> Reach gate on the projection distance, falling from `reach_on` to `reach_off`
      type(moist_math_smoothstep_sigmoid_type) :: f_reach
      type(branch_weight_type) :: branch_weight
      !> Seeds per anchor in the first round of the extended search, doubled per round
      integer :: nearest_count = default_nearest_count
      !> Keep searching after the first accepted point of an anchor without a local
      !> point; off, the nearest point of that round is the anchor's only row
      logical :: confirm_nearest = .true.
      !> Pass-1 (local), extended-search and pass-2 (seeded) batches
      type(projection_batch_type) :: pass1, passn, pass2
      !> Anchors in the extended search: pass-1 run, first, last and next seed
      integer, allocatable :: near_run(:), near_first(:), near_last(:), near_next(:)
      !> Their accepted local column in pass 1 (0: none) and latest run in `passn` (0: none)
      integer, allocatable :: near_jloc(:), near_tail(:)
      !> Whether the anchor is seeded by `fallback`
      logical, allocatable :: near_fallback(:)
      !> Seeder reach of each such anchor's seeds, the radius its points ask for
      !> and the radius of the ball its seeds cover (Bohr)
      real(wp), allocatable :: near_reach(:), near_radius(:), near_ball(:)
      !> Nearest-first seeds of those anchors (3, :)
      real(wp), allocatable :: near_seeds(:, :)
      !> Scratches
      type(newton_scratch_type) :: nscratch
      type(finalize_scratch_type) :: fscratch
      !> Rows appended by every chunk this worker processed
      type(projection_rows_type) :: rows
      !> Anchors seeded, and seeds proposed, since `init`
      integer :: n_seeded = 0
      integer :: n_seeds = 0
      !> Anchors without an accepted local point, seeded since `init`
      integer :: n_escalated = 0
      !> Seeded columns stopped as duplicates of an accepted point, since `init`
      integer :: n_captured = 0
      !> Anchors searched beyond the window of a fully open local point, since `init`
      integer :: n_extended = 0
      !> Largest search radius of an anchor since `init` (Bohr)
      real(wp) :: pdist_search_max = 0.0_wp
      !> Anchors searched to the reach from the start since `init`
      integer :: n_full = 0
      !> Anchors whose seeds all failed and were solved again, since `init`
      integer :: n_restarted = 0
   contains
      !> Configure the worker from the cavity parameters
      procedure :: init => worker_init
      !> Run the pipeline on one chunk of anchors
      procedure :: process_chunk => worker_process_chunk
      !> Single-anchor entry: project one anchor and return its rows
      procedure :: project_anchor => worker_project_anchor
   end type projection_worker_type

contains

   !> Configure the worker from the cavity parameters
   !>
   !> - the LSF clone may borrow read-only data of `lsf_model`, which must outlive
   !>   the worker unchanged
   !>
   !> @param[inout] self      Worker
   !> @param[in]    param     DROP parameters (validated)
   !> @param[in]    lsf_model Level set, already updated for the geometry; cloned
   !> @param[in]    centers   Atom centers (3, nat)
   !> @param[out]   error     Invalid seeder configuration
   !> @param[in]    debug     Per-anchor seeder traces
   !> @param[in]    unit      Output unit of the seeder traces, standard output when absent
   subroutine worker_init(self, param, lsf_model, centers, error, debug, unit)
      class(projection_worker_type), intent(inout) :: self
      type(moist_cavity_drop_parameters_type), intent(in) :: param
      class(moist_cavity_drop_lsf_type), intent(in), target :: lsf_model
      real(wp), intent(in) :: centers(:, :)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in), optional :: debug
      integer, intent(in), optional :: unit

      if (.not. is_valid_proj_level(param%proj_level)) then
         call fatal_error(error, "DROP projection: level must be 1, 3 or 4")
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
      self%solve_restart = self%solve_seeded
      self%solve_restart%max_cheap_steps = 0
      self%polish = param%polish
      self%finalize = param%finalize
      self%phi_alpha = param%phi_alpha

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

      if (allocated(self%seeder)) deallocate (self%seeder)
      if (allocated(self%fallback)) deallocate (self%fallback)
      select case (self%level)
      case (1)
         call new_octree_seeder(self%fallback, param%seed, error, debug, unit)
      case (3, 4)
         call new_octree_seeder(self%seeder, param%seed, error, debug, unit)
      case default
         error stop "DROP projection worker: level other than 1, 3, 4 passed the level check"
      end select
      if (allocated(error)) return

      call self%pass1%init(param%batch_size, error)
      if (allocated(error)) return
      call self%passn%init(param%batch_size, error)
      if (allocated(error)) return
      call self%pass2%init(param%batch_size, error)
      if (allocated(error)) return
      call self%rows%clear()
      self%n_seeded = 0
      self%n_seeds = 0
      self%n_escalated = 0
      self%n_captured = 0
      self%n_extended = 0
      self%pdist_search_max = 0.0_wp
      self%n_full = 0
      self%n_restarted = 0
   end subroutine worker_init

   !> Allocate and configure an octree seeder
   !>
   !> @param[out] seeder        Octree seeder slot
   !> @param[in]  seed          Seed-step settings of the level
   !> @param[out] error         Invalid octree configuration
   !> @param[in]  debug         Per-anchor search report
   !> @param[in]  unit          Output unit of the report
   subroutine new_octree_seeder(seeder, seed, error, debug, unit)
      class(projection_seeder_type), allocatable, intent(out) :: seeder
      type(drop_seed_parameters_type), intent(in) :: seed
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in), optional :: debug
      integer, intent(in), optional :: unit

      allocate (octree_seeder_type :: seeder)
      select type (seeder)
      type is (octree_seeder_type)
         if (present(debug)) seeder%debug = debug
         if (present(unit)) seeder%unit = unit
         call seeder%init(error, seed_size=seed%seed_size, max_boxes=seed%max_boxes, &
                          max_survivors=seed%max_survivors, max_depth=seed%max_depth, &
                          refine_levels=seed%refine_levels)
      end select
   end subroutine new_octree_seeder

   !> Run the pipeline on one chunk of anchors
   !>
   !> @param[inout] self    Worker
   !> @param[in]    ctx     Run context, read only (polish traces)
   !> @param[in]    slots   Anchor slots, written into the rows (n)
   !> @param[in]    anchors Anchor positions (3, n)
   !> @param[in]    owners  Owning atoms (n)
   !> @param[in]    wlebs   Anchor quadrature weights (n)
   !> @param[in]    fs      Anchor switching values (n); below `wleb_cut` the anchor is skipped
   !> @param[out]   nrows   Rows appended per anchor (n)
   !> @param[out]   lost    Anchor kept no accepted point (n)
   !> @param[out]   error   Build-fatal failure
   !> @param[out]   unseeded Lost anchor whose searches up to the reach gave no seed (n)
   subroutine worker_process_chunk(self, ctx, slots, anchors, owners, wlebs, fs, nrows, lost, error, unseeded)
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

      integer :: k, r, j, jacc, icol, nnear, nnear_seeds
      real(wp) :: gate

      nrows = 0
      lost = .false.
      if (present(unseeded)) unseeded = .false.

      !> Pass 1: local Newton from every anchor above the switching cut
      call self%pass1%clear()
      do k = 1, size(slots)
         if (fs(k) < self%finalize%wleb_cut) cycle
         call self%pass1%begin_run(k, anchors(:, k), owners(k))
         call self%pass1%append(anchors(:, k), 0.0_wp, icol)
      end do
      if (self%pass1%n == 0) return
      call newton_solve_batch(self%pass1, self%lsf, self%solve, self%phi_alpha, self%nscratch, error)
      if (allocated(error)) return
      call polish_batch(ctx, self%pass1, self%lsf, self%polish, self%solve, self%phi_alpha, &
                        self%nscratch, error)
      if (allocated(error)) return

      ! Seeding; only a fully open local point bounds the search by its window
      call self%pass2%clear()
      call self%passn%clear()
      nnear = 0
      nnear_seeds = 0
      do r = 1, self%pass1%nrun
         k = self%pass1%run_anchor(r)
         j = self%pass1%run_start(r)
         jacc = 0
         if (self%pass1%status(j) == point_accepted) jacc = j
         if (jacc == 0) then
            self%n_escalated = self%n_escalated + 1
            call queue_search(allocated(self%fallback), r, k, 0)
         else if (.not. allocated(self%seeder)) then
            call finish_anchor(self%pass1, j, j, k)
         else
            call column_gate(self%pass1, jacc, k, gate)
            if (allocated(error)) return
            if (gate == 1.0_wp .and. .not. self%full_reach) then
               call seed_anchor(self%seeder, jacc, r, k)
            else
               call queue_search(.false., r, k, jacc)
            end if
         end if
         if (allocated(error)) return
      end do
      call solve_queue()
      if (allocated(error)) return
      call flush_pass2()

   contains

      !> Request of the anchor of pass-1 run `r1`: the whole ball of radius `radius`
      !>
      !> @param[in]  r1      Pass-1 run index
      !> @param[in]  k1      Anchor index within the chunk
      !> @param[in]  radius  Search radius (Bohr)
      !> @param[out] request Request
      subroutine anchor_request(r1, k1, radius, request)
         integer, intent(in) :: r1, k1
         real(wp), intent(in) :: radius
         type(seed_request_type), intent(out) :: request

         request%index = slots(k1)
         request%anchor = anchors(:, k1)
         request%owner_center = self%centers(:, owners(k1))
         call self%pass1%anchor_jet(r1, request%lsf0_anchor, request%grad_anchor)
         request%pdist2_slack = self%branch_weight%x_max**2
         request%pdist_slack = self%branch_weight%x_max
         request%pdist_soft = self%branch_weight%soft
         request%pdist_fallback = self%reach_off
         request%pdist_start = radius
         request%fixed_radius = .true.
      end subroutine anchor_request

      !> Search radius behind a fully open point at distance `pdist`
      !>
      !> @param[in] pdist  Distance of the point (Bohr)
      !> @return    radius Largest distance inside its branch window, at most the reach (Bohr)
      pure function window_radius(pdist) result(radius)
         real(wp), intent(in) :: pdist
         real(wp) :: radius

         radius = min(self%reach_off, self%branch_weight%pdist_max(pdist))
      end function window_radius

      !> Finalize record of an anchor of the chunk
      !>
      !> @param[in]  k2  Anchor index within the chunk
      !> @param[out] rec Per-anchor scalars
      subroutine anchor_record(k2, rec)
         integer, intent(in) :: k2
         type(finalize_anchor_type), intent(out) :: rec

         rec%index = slots(k2)
         rec%xyz = anchors(:, k2)
         rec%owner = owners(k2)
         rec%owner_center = self%centers(:, owners(k2))
         rec%wleb = wlebs(k2)
         rec%f = fs(k2)
      end subroutine anchor_record

      !> Gate of an accepted column, the one the finalize will see
      !>
      !> @param[in]  batch  Batch holding the column
      !> @param[in]  jcol   Accepted column
      !> @param[in]  k2     Anchor index within the chunk
      !> @param[out] gate_j Gate
      subroutine column_gate(batch, jcol, k2, gate_j)
         type(projection_batch_type), intent(in) :: batch
         integer, intent(in) :: jcol, k2
         real(wp), intent(out) :: gate_j
         type(finalize_anchor_type) :: rec
         real(wp) :: jac

         call anchor_record(k2, rec)
         call point_gate(self%finalize, self%phi_alpha, batch, jcol, rec, self%f_crit, self%f_foc, self%f_hor, &
                         self%f_wleb, self%f_reach, jac, gate_j, error)
      end subroutine column_gate

      !> Seed the window ball of a fully open local point and open its pass-2 run,
      !> flushing the batch first when full
      !>
      !> @param[inout] seeder Seeder
      !> @param[in]    jacc   Pass-1 column of the accepted local point
      !> @param[in]    r1     Pass-1 run index
      !> @param[in]    k1     Anchor index within the chunk
      subroutine seed_anchor(seeder, jacc, r1, k1)
         class(projection_seeder_type), intent(inout) :: seeder
         integer, intent(in) :: jacc, r1, k1
         type(seed_request_type) :: request
         real(wp) :: radius
         integer :: n_new, i

         radius = window_radius(self%pass1%pdist(jacc))
         self%pdist_search_max = max(self%pdist_search_max, radius)
         call anchor_request(r1, k1, radius, request)
         call seeder%seed(self%lsf, request, error)
         if (allocated(error)) return
         self%n_seeded = self%n_seeded + 1
         self%n_seeds = self%n_seeds + seeder%n_points

         n_new = seeder%n_points + 1
         if (self%pass2%would_overflow(n_new)) then
            call flush_pass2()
            if (allocated(error)) return
         end if
         call self%pass2%begin_run(k1, anchors(:, k1), owners(k1))
         call self%pass2%append_copy(self%pass1, jacc, icol)
         ! Seeds start with a zero multiplier
         do i = 1, seeder%n_points
            call self%pass2%append(seeder%points(:, i), 0.0_wp, icol)
         end do
      end subroutine seed_anchor

      !> Queue the anchor of pass-1 run `r1` for the extended search
      !>
      !> - first ball: the window ball of the local point, else of the anchor's
      !>   Newton length `|S|/|grad S|`, a guess the search outgrows when wrong
      !>
      !> @param[in] use_fallback Seed with the fallback seeder
      !> @param[in] r1           Pass-1 run index
      !> @param[in] k1           Anchor index within the chunk
      !> @param[in] jacc         Pass-1 column of the accepted local point, 0 for none
      subroutine queue_search(use_fallback, r1, k1, jacc)
         logical, intent(in) :: use_fallback
         integer, intent(in) :: r1, k1, jacc
         real(wp) :: s_anchor, g_anchor(3), ball

         if (jacc > 0) then
            ball = window_radius(self%pass1%pdist(jacc))
         else
            call self%pass1%anchor_jet(r1, s_anchor, g_anchor)
            ball = self%reach_off
            if (norm2(g_anchor) > 0.0_wp) ball = window_radius(min(self%reach_off, abs(s_anchor)/norm2(g_anchor)))
         end if
         if (self%full_reach) ball = self%reach_off
         ball = max(ball, min_ball_fraction*self%reach_off)
         if (ball >= self%reach_off) self%n_full = self%n_full + 1

         nnear = nnear + 1
         call reserve_queue(nnear, nnear_seeds)
         self%near_run(nnear) = r1
         self%near_jloc(nnear) = jacc
         self%near_tail(nnear) = 0
         self%near_fallback(nnear) = use_fallback
         self%near_radius(nnear) = self%reach_off
         self%near_ball(nnear) = ball
         self%n_seeded = self%n_seeded + 1
         self%n_extended = self%n_extended + 1
         call seed_ball(nnear, k1)
      end subroutine queue_search

      !> Seed the ball of a queued anchor, growing it while it holds no seed
      !>
      !> - the seeds replace the anchor's unsolved ones
      !>
      !> @param[in] e  Queue index of the anchor
      !> @param[in] k1 Anchor index within the chunk
      subroutine seed_ball(e, k1)
         integer, intent(in) :: e, k1
         type(seed_request_type) :: request
         integer :: n

         ! TODO: a grown ball is seeded whole, so its inner part is probed and
         ! solved again (12 % of the seeds on c8 at level 3); seed only the
         ! shell, skipping seeds within `r_old - h_new` of the anchor, and keep
         ! the earlier seeds ahead of it so `restart_run` still sees all of them
         do
            call anchor_request(self%near_run(e), k1, self%near_ball(e), request)
            request%nearest_first = .true.
            if (self%near_fallback(e)) then
               call self%fallback%seed(self%lsf, request, error)
               if (allocated(error)) return
               n = self%fallback%n_points
               call reserve_queue(nnear, nnear_seeds + n)
               ! `points` stays unallocated until the first seed
               if (n > 0) self%near_seeds(:, nnear_seeds + 1:nnear_seeds + n) = self%fallback%points(:, 1:n)
               self%near_reach(e) = self%fallback%reach
            else
               call self%seeder%seed(self%lsf, request, error)
               if (allocated(error)) return
               n = self%seeder%n_points
               call reserve_queue(nnear, nnear_seeds + n)
               ! `points` stays unallocated until the first seed
               if (n > 0) self%near_seeds(:, nnear_seeds + 1:nnear_seeds + n) = self%seeder%points(:, 1:n)
               self%near_reach(e) = self%seeder%reach
            end if
            self%near_first(e) = nnear_seeds + 1
            self%near_next(e) = nnear_seeds + 1
            nnear_seeds = nnear_seeds + n
            self%near_last(e) = nnear_seeds
            if (n > 0 .or. self%near_ball(e) >= self%near_radius(e)) exit
            self%near_ball(e) = min(self%near_radius(e), 2.0_wp*self%near_ball(e))
         end do
      end subroutine seed_ball

      !> Grow the search queue to `n` anchors and `m` seeds, keeping its contents
      !>
      !> @param[in] n Anchors to hold
      !> @param[in] m Seeds to hold
      subroutine reserve_queue(n, m)
         integer, intent(in) :: n, m
         real(wp), allocatable :: tmp_seeds(:, :)
         logical, allocatable :: tmp_flag(:)
         integer :: cap

         if (.not. allocated(self%near_run)) then
            allocate (self%near_run(64), self%near_last(64), self%near_next(64), self%near_reach(64))
            allocate (self%near_first(64))
            allocate (self%near_jloc(64), self%near_tail(64), self%near_radius(64), self%near_ball(64))
            allocate (self%near_fallback(64))
            allocate (self%near_seeds(3, 1024))
         end if
         if (n > size(self%near_run)) then
            cap = max(n, 2*size(self%near_run))
            call grow_real(self%near_reach, cap, n - 1)
            call grow_real(self%near_radius, cap, n - 1)
            call grow_real(self%near_ball, cap, n - 1)
            call grow_int(self%near_run, cap, n - 1)
            call grow_int(self%near_first, cap, n - 1)
            call grow_int(self%near_last, cap, n - 1)
            call grow_int(self%near_next, cap, n - 1)
            call grow_int(self%near_jloc, cap, n - 1)
            call grow_int(self%near_tail, cap, n - 1)
            allocate (tmp_flag(cap))
            tmp_flag(1:n - 1) = self%near_fallback(1:n - 1)
            call move_alloc(tmp_flag, self%near_fallback)
         end if
         if (m > size(self%near_seeds, 2)) then
            cap = max(m, 2*size(self%near_seeds, 2))
            allocate (tmp_seeds(3, cap))
            tmp_seeds(:, 1:nnear_seeds) = self%near_seeds(:, 1:nnear_seeds)
            call move_alloc(tmp_seeds, self%near_seeds)
         end if
      end subroutine reserve_queue

      !> Extended search of every queued anchor, in rounds batched across anchors
      !>
      !> - round `i` solves the next `nearest_count 2^(i-1)` seeds of every
      !>   pending anchor that lie within the seeder's reach of its radius
      !> - after the round the radius of an anchor shrinks to the window behind
      !>   its nearest fully open accepted point
      !> - an anchor with no unsolved seed left inside gets a larger ball while
      !>   its ball is short of its radius, and else is finalized from its
      !>   latest run, which holds every point accepted for it
      subroutine solve_queue()
         integer :: take, npend, nkeep, ip, e, k1, rn, s, s_last, jn, j_first, j_last, col0, run0, jbest
         integer, allocatable :: pend(:)
         real(wp) :: limit, gate_j
         logical :: more

         if (nnear == 0) return
         ! A round without seeds would leave every waiting anchor waiting forever
         if (self%nearest_count < 1) then
            call fatal_error(error, "DROP projection: nearest_count must be positive")
            return
         end if
         allocate (pend(nnear))
         do e = 1, nnear
            pend(e) = e
         end do
         npend = nnear
         take = self%nearest_count
         do while (npend > 0)
            col0 = self%passn%n + 1
            run0 = self%passn%nrun + 1
            do ip = 1, npend
               e = pend(ip)
               k1 = self%pass1%run_anchor(self%near_run(e))
               limit = self%near_radius(e) + self%near_reach(e)
               s_last = self%near_next(e) - 1
               do s = self%near_next(e), min(self%near_last(e), self%near_next(e) + take - 1)
                  if (norm2(self%near_seeds(:, s) - anchors(:, k1)) > limit) exit
                  s_last = s
               end do
               if (s_last < self%near_next(e)) cycle

               ! Points accepted so far lead the run: kept, and captors of the new seeds
               ! TODO: superseded runs stay in passn and accepted points are recopied every round;
               ! growth ~ rounds x accepted per anchor, undocumented and never reclaimed
               j_first = 1
               j_last = 0
               if (self%near_tail(e) > 0) then
                  j_first = self%passn%run_start(self%near_tail(e))
                  j_last = j_first + self%passn%run_len(self%near_tail(e)) - 1
               end if
               call self%passn%begin_run(e, anchors(:, k1), owners(k1))
               if (self%near_tail(e) > 0) then
                  do jn = j_first, j_last
                     if (self%passn%status(jn) == point_accepted) call self%passn%append_own(jn, icol)
                  end do
               else if (self%near_jloc(e) > 0) then
                  call self%passn%append_copy(self%pass1, self%near_jloc(e), icol)
               end if
               self%near_tail(e) = self%passn%nrun
               do s = self%near_next(e), s_last
                  call self%passn%append(self%near_seeds(:, s), 0.0_wp, icol)
               end do
               self%n_seeds = self%n_seeds + s_last - self%near_next(e) + 1
               self%near_next(e) = s_last + 1
            end do

            if (self%passn%n >= col0) then
               call newton_solve_batch(self%passn, self%lsf, self%solve_seeded, self%phi_alpha, self%nscratch, &
                                       error, first=col0)
               if (allocated(error)) return
               call polish_batch(ctx, self%passn, self%lsf, self%polish, self%solve_seeded, self%phi_alpha, &
                                 self%nscratch, error, first=col0)
               if (allocated(error)) return
               self%n_captured = self%n_captured + count(self%passn%status(col0:self%passn%n) == point_captured)
            end if

            ! Radius behind the nearest fully open point of each anchor
            if (.not. self%full_reach) then
               do rn = run0, self%passn%nrun
                  e = self%passn%run_anchor(rn)
                  k1 = self%pass1%run_anchor(self%near_run(e))
                  do jn = self%passn%run_start(rn), self%passn%run_start(rn) + self%passn%run_len(rn) - 1
                     if (self%passn%status(jn) /= point_accepted) cycle
                     call column_gate(self%passn, jn, k1, gate_j)
                     if (allocated(error)) return
                     if (gate_j == 1.0_wp) then
                        self%near_radius(e) = min(self%near_radius(e), window_radius(self%passn%pdist(jn)))
                     end if
                  end do
               end do
            end if

            nkeep = 0
            do ip = 1, npend
               e = pend(ip)
               k1 = self%pass1%run_anchor(self%near_run(e))
               more = self%near_next(e) <= self%near_last(e)
               if (more) then
                  more = norm2(self%near_seeds(:, self%near_next(e)) - anchors(:, k1)) &
                         <= self%near_radius(e) + self%near_reach(e)
               end if
               ! Without confirmation the nearest point of the first successful round stands alone
               jbest = 0
               if (.not. self%confirm_nearest .and. self%near_jloc(e) == 0 .and. self%near_tail(e) > 0) then
                  rn = self%near_tail(e)
                  do jn = self%passn%run_start(rn), self%passn%run_start(rn) + self%passn%run_len(rn) - 1
                     if (self%passn%status(jn) /= point_accepted) cycle
                     if (jbest == 0) then
                        jbest = jn
                     else if (self%passn%pdist(jn) < self%passn%pdist(jbest)) then
                        jbest = jn
                     end if
                  end do
               end if
               if (jbest > 0) then
                  call finish_anchor(self%passn, jbest, jbest, k1)
               else if (more) then
                  nkeep = nkeep + 1
                  pend(nkeep) = e
               else if (self%near_ball(e) < self%near_radius(e)) then
                  ! To the radius a fully open point set, else twice as far
                  if (self%near_radius(e) < self%reach_off) then
                     self%near_ball(e) = self%near_radius(e)
                  else
                     self%near_ball(e) = min(self%near_radius(e), 2.0_wp*self%near_ball(e))
                  end if
                  call seed_ball(e, k1)
                  nkeep = nkeep + 1
                  pend(nkeep) = e
               else
                  call finish_queued(e, k1)
               end if
               if (allocated(error)) return
            end do
            npend = nkeep
            take = 2*take
         end do
      end subroutine solve_queue

      !> Finalize a queued anchor from its latest run, or from pass 1 without one
      !>
      !> @param[in] e  Queue index of the anchor
      !> @param[in] k1 Anchor index within the chunk
      subroutine finish_queued(e, k1)
         integer, intent(in) :: e, k1
         integer :: rn, j1

         self%pdist_search_max = max(self%pdist_search_max, self%near_ball(e))
         if (self%near_tail(e) > 0) then
            rn = self%near_tail(e)
            j1 = self%passn%run_start(rn)
            if (.not. any(self%passn%status(j1:j1 + self%passn%run_len(rn) - 1) == point_accepted)) then
               call restart_run(e, k1)
               if (allocated(error)) return
               rn = self%near_tail(e)
            end if
            call finish_anchor(self%passn, self%passn%run_start(rn), &
                               self%passn%run_start(rn) + self%passn%run_len(rn) - 1, k1)
         else
            j1 = self%pass1%run_start(self%near_run(e))
            call finish_anchor(self%pass1, j1, j1, k1)
            ! No seed was ever solved: no sign change within the reach
            if (present(unseeded)) unseeded(k1) = lost(k1)
         end if
      end subroutine finish_queued

      !> Solve the seeds of a queued anchor again, with the full system from the start
      !>
      !> - for an anchor whose seeds all failed: the cheap step puts a seed at
      !>   `a + c grad S`, which leaves a thin neck around an anchor on its axis
      !> - the seeds of its last ball, nearest first, form its new latest run
      !>
      !> @param[in] e  Queue index of the anchor
      !> @param[in] k1 Anchor index within the chunk
      subroutine restart_run(e, k1)
         integer, intent(in) :: e, k1
         integer :: s, col0
         real(wp) :: limit

         limit = self%near_radius(e) + self%near_reach(e)
         col0 = self%passn%n + 1
         call self%passn%begin_run(e, anchors(:, k1), owners(k1))
         self%near_tail(e) = self%passn%nrun
         do s = self%near_first(e), self%near_last(e)
            if (norm2(self%near_seeds(:, s) - anchors(:, k1)) > limit) exit
            call self%passn%append(self%near_seeds(:, s), 0.0_wp, icol)
         end do
         self%n_restarted = self%n_restarted + 1
         if (self%passn%n < col0) return
         call newton_solve_batch(self%passn, self%lsf, self%solve_restart, self%phi_alpha, self%nscratch, &
                                 error, first=col0)
         if (allocated(error)) return
         call polish_batch(ctx, self%passn, self%lsf, self%polish, self%solve_restart, self%phi_alpha, &
                           self%nscratch, error, first=col0)
         if (allocated(error)) return
         self%n_captured = self%n_captured + count(self%passn%status(col0:self%passn%n) == point_captured)
      end subroutine restart_run

      !> Solve, polish and finalize every run of the pass-2 batch
      subroutine flush_pass2()
         integer :: r2, k2

         if (self%pass2%nrun == 0) return
         if (self%pass2%n > 0) then
            call newton_solve_batch(self%pass2, self%lsf, self%solve_seeded, self%phi_alpha, self%nscratch, error)
            if (allocated(error)) return
            call polish_batch(ctx, self%pass2, self%lsf, self%polish, self%solve_seeded, self%phi_alpha, &
                              self%nscratch, error)
            if (allocated(error)) return
            self%n_captured = self%n_captured + count(self%pass2%status(1:self%pass2%n) == point_captured)
         end if
         do r2 = 1, self%pass2%nrun
            k2 = self%pass2%run_anchor(r2)
            call finish_anchor(self%pass2, self%pass2%run_start(r2), &
                               self%pass2%run_start(r2) + self%pass2%run_len(r2) - 1, k2)
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

         call anchor_record(k2, rec)
         call finalize_anchor(self%finalize, self%phi_alpha, batch, first, last, rec, self%f_crit, &
                              self%f_foc, self%f_hor, self%f_wleb, self%f_reach, self%branch_weight, self%rows, &
                              npoint, nrows(k2), self%fscratch, error)
         lost(k2) = npoint == 0
      end subroutine finish_anchor

   end subroutine worker_process_chunk

   !> Project one anchor and return its rows (the single-anchor entry for tests)
   !>
   !> - the worker's row buffer is cleared first, so the outputs are exactly
   !>   this anchor's surviving branches in canonical order
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
   !> @param[out]   n_seeds   Seeds the seeder proposed for this anchor
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
      if (present(points)) points = self%rows%xyz(:, 1:n_points)
      if (present(normals)) normals = self%rows%normal(:, 1:n_points)
      if (present(lambda)) lambda = self%rows%lambda(1:n_points)
      if (present(pdist)) pdist = self%rows%pdist(1:n_points)
      if (present(phi)) phi = self%rows%phi(1:n_points)
      if (present(weights)) weights = self%rows%wbranch(1:n_points)
   end subroutine worker_project_anchor

   !> Grow an integer vector to `cap`, keeping its first `keep` entries
   !>
   !> @param[inout] a    Vector
   !> @param[in]    cap  New length
   !> @param[in]    keep Leading entries to keep
   subroutine grow_int(a, cap, keep)
      integer, allocatable, intent(inout) :: a(:)
      integer, intent(in) :: cap, keep
      integer, allocatable :: tmp(:)

      allocate (tmp(cap))
      tmp(1:keep) = a(1:keep)
      call move_alloc(tmp, a)
   end subroutine grow_int

   !> Grow a real vector to `cap`, keeping its first `keep` entries
   !>
   !> @param[inout] a    Vector
   !> @param[in]    cap  New length
   !> @param[in]    keep Leading entries to keep
   subroutine grow_real(a, cap, keep)
      real(wp), allocatable, intent(inout) :: a(:)
      integer, intent(in) :: cap, keep
      real(wp), allocatable :: tmp(:)

      allocate (tmp(cap))
      tmp(1:keep) = a(1:keep)
      call move_alloc(tmp, a)
   end subroutine grow_real

end module moist_cavity_drop_projection_worker
