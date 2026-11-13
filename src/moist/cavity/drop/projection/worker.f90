!> Per-thread worker of the batched DROP projection
!>
!> One worker runs the whole per-anchor pipeline on a chunk of anchors:
!>
!>   1. local Newton from every anchor (pass-1 batch), optimality polish
!>   2. global search: every anchor at levels 3, 4, none at level 1
!>   3. seeder on the searched anchors; their seeds, plus the accepted local
!>      foot, form one run each in the pass-2 batch, flushed when full
!>   4. nearest-foot search on anchors without an accepted local foot (below)
!>   5. batched Newton and polish on the seeds; a seed whose Newton prediction
!>      lands on an accepted foot of its anchor stops as a duplicate
!>   6. per-anchor finalize (dedupe, canonical order, CP Jacobian, switching,
!>      branch weights, cut) appending rows to the worker's row buffer
!>
!> Nearest-foot search of an anchor without an accepted local foot:
!>
!>   - the fallback seeder (else the level's seeder) proposes its seeds
!>     nearest first; the nearest `nearest_count` of every such anchor of the
!>     chunk are solved in one batch, and an anchor with none accepted takes the
!>     next ones, twice as many per round, until a foot is accepted or its seeds
!>     run out
!>   - the accepted foot closest to the anchor then stands in for the local
!>     foot: its pass-2 run takes the round's accepted feet and every unsolved
!>     seed within the seeder's `reach` of `sqrt(rho^2 + slack)`, the seeds a
!>     search from that foot would propose, so a closer foot and every branch
!>     in the window are still found without searching again
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
!> - an anchor without an accepted local foot goes to the `fallback` seeder
!>   (level 1: level-3 octree settings), else to the level's seeder; with
!>   nothing accepted at the end it is lost, reported through `lost`, and
!>   appends no row
!> - a surviving leaf without sign-change evidence is dropped by the coarse
!>   octree and refined once by the fine one
!> - `nearest_count` is a starting guess, not a correctness parameter: the
!>   second search finds a closer foot the first round missed
!> - LSF failures, an exhausted octree budget, an anchor on its owner center
!>   and a singular tangent block are build-fatal and come back as `error`
module moist_cavity_drop_projection_worker
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type, is_valid_proj_level
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_projection_batch, only: projection_batch_type, point_accepted, point_captured
   use moist_cavity_drop_projection_newton, only: newton_config_type, newton_scratch_type, &
                                                 newton_solve_batch
   use moist_cavity_drop_projection_polish, only: polish_config_type, polish_batch
   use moist_cavity_drop_projection_finalize, only: finalize_config_type, finalize_anchor_type, &
                                                   finalize_scratch_type, projection_rows_type, &
                                                   finalize_anchor
   use moist_cavity_drop_projection_seed, only: projection_seeder_type, seed_request_type
   use moist_cavity_drop_projection_seed_octree, only: octree_seeder_type
   implicit none(type, external)
   private

   public :: projection_worker_type

   !> Octree search radius around an anchor when neither a local foot nor the
   !> surface trace bounds the search (Bohr)
   real(wp), parameter :: fallback_search_radius = 15.0_wp
   !> Seeds per anchor in the first round of the nearest-foot search
   integer, parameter :: default_nearest_count = 4
   !> Refinement levels of an evidence-free leaf in the coarse octree (level 3, fallback of 1)
   integer, parameter :: coarse_refine_levels = 0
   !> Refinement levels of an evidence-free leaf in the fine octree (level 4)
   integer, parameter :: fine_refine_levels = 1

   !> Per-thread projection worker
   type :: projection_worker_type
      !> Projection level (1, 3 or 4)
      integer :: level = 3
      !> Thread-local level set clone
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      !> Atom centers (3, nat)
      real(wp), allocatable :: centers(:, :)
      !> Stage settings
      type(newton_config_type) :: ncfg
      type(polish_config_type) :: pcfg
      type(finalize_config_type) :: fcfg
      !> Newton settings of the seeded pass: `ncfg` with duplicate capture
      type(newton_config_type) :: ncfg_seeded
      !> Seeder of the level, searching every anchor; unallocated at level 1
      class(projection_seeder_type), allocatable :: seeder
      !> Seeder of anchors without an accepted local foot; unallocated: `seeder`
      class(projection_seeder_type), allocatable :: fallback
      !> Switching functions and branch weight model
      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc, f_wleb
      type(branch_weight_type) :: branch_weight
      !> Seeds per anchor in the first nearest-foot round, doubled per round
      integer :: nearest_count = default_nearest_count
      !> Search again around the nearest foot; off, that foot is the anchor's only row
      logical :: confirm_nearest = .true.
      !> Pass-1 (local), nearest-foot and pass-2 (seeded) batches
      type(projection_batch_type) :: pass1, passn, pass2
      !> Anchors awaiting the nearest-foot search: pass-1 run, last and next seed
      integer, allocatable :: near_run(:), near_last(:), near_next(:)
      !> Seeder reach of each waiting anchor's seeds (Bohr)
      real(wp), allocatable :: near_reach(:)
      !> Nearest-first seeds of the waiting anchors (3, :)
      real(wp), allocatable :: near_seeds(:, :)
      !> Scratches
      type(newton_scratch_type) :: nscratch
      type(finalize_scratch_type) :: fscratch
      !> Rows appended by every chunk this worker processed
      type(projection_rows_type) :: rows
      !> Anchors seeded, and seeds proposed, since `init`
      integer :: n_seeded = 0
      integer :: n_seeds = 0
      !> Anchors without an accepted local foot, seeded since `init`
      integer :: n_escalated = 0
      !> Seeded columns stopped as duplicates of an accepted foot, since `init`
      integer :: n_captured = 0
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
   !> @param[inout] self      Worker
   !> @param[in]    param     DROP parameters (validated)
   !> @param[in]    lsf_model Level set, already updated for the geometry; cloned
   !> @param[in]    centers   Atom centers (3, nat)
   !> @param[out]   error     Invalid seeder configuration
   !> @param[in]    verbosity Run-context verbosity
   !> @param[in]    debug     Per-anchor debug traces
   !> @param[in]    unit      Output unit of the traces, standard output when absent
   subroutine worker_init(self, param, lsf_model, centers, error, verbosity, debug, unit)
      class(projection_worker_type), intent(inout) :: self
      type(moist_cavity_drop_parameters_type), intent(in) :: param
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf_model
      real(wp), intent(in) :: centers(:, :)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: verbosity
      logical, intent(in), optional :: debug
      integer, intent(in), optional :: unit

      !> Octree leaf edge of the level
      real(wp) :: seed_size
      !> Octree refinement levels of the level
      integer :: refine_levels

      if (.not. is_valid_proj_level(param%proj_level)) then
         call fatal_error(error, "DROP projection: level must be 1, 3 or 4")
         return
      end if
      self%level = param%proj_level

      if (allocated(self%lsf)) deallocate (self%lsf)
      allocate (self%lsf, source=lsf_model)
      call self%lsf%set_max_deriv(2)
      self%centers = centers

      self%ncfg%phi_alpha = param%phi_alpha
      self%ncfg%tol = param%proj_tol
      self%ncfg%max_iter = param%proj_maxiter
      self%ncfg%switch_fnorm = param%newton_switch_fnorm
      self%ncfg%max_cheap_steps = param%newton_max_cheap_steps
      self%ncfg%max_step = param%newton_max_step
      self%ncfg%merit_margin = param%newton_merit_margin
      self%ncfg_seeded = self%ncfg
      self%ncfg_seeded%capture_duplicates = .true.

      self%pcfg%phi_alpha = param%phi_alpha
      if (present(verbosity)) self%pcfg%verbosity = verbosity
      if (present(debug)) self%pcfg%debug = debug
      if (present(unit)) self%pcfg%unit = unit

      self%fcfg%phi_alpha = param%phi_alpha
      self%fcfg%branch_sep_cut = param%branch_sep_cut
      self%fcfg%rho2_slack = param%branch_rho2_slack()
      self%fcfg%wleb_cut = param%wleb_cut
      self%fcfg%wleb_prune_level = param%wleb_prune_level

      call new_smoothstep_sigmoid(self%f_crit, param%w_0ls_from, param%w_0ls_to, &
                                  p_hi=param%w_0ls_p, a_hi=param%w_0ls_a, p_lo=param%w_0ls_p, &
                                  a_lo=param%w_0ls_a)
      call new_smoothstep_sigmoid(self%f_foc, param%w_0tra_from, param%w_0tra_to, &
                                  p_hi=param%w_0ls_p, a_hi=param%w_0ls_a, p_lo=param%w_0ls_p, &
                                  a_lo=param%w_0ls_a)
      if (param%wleb_prune_level > 0) then
         call new_smoothstep_sigmoid(self%f_wleb, param%wleb_prune_from, param%wleb_prune_to)
      else
         call new_smoothstep_sigmoid(self%f_wleb, 0.0_wp, 1.0_wp)
      end if
      call self%branch_weight%init(param%branch_weight_s)

      if (allocated(self%seeder)) deallocate (self%seeder)
      if (allocated(self%fallback)) deallocate (self%fallback)
      select case (self%level)
      case (1)
         call new_octree_seeder(self%fallback, param, param%octree_coarse_seed_size, &
                                coarse_refine_levels, error, debug, unit)
      case (3, 4)
         seed_size = param%octree_seed_size
         refine_levels = fine_refine_levels
         if (self%level == 3) then
            seed_size = param%octree_coarse_seed_size
            refine_levels = coarse_refine_levels
         end if
         call new_octree_seeder(self%seeder, param, seed_size, refine_levels, error, debug, unit)
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
   end subroutine worker_init

   !> Allocate and configure an octree seeder
   !>
   !> @param[out] seeder        Octree seeder slot
   !> @param[in]  param         DROP parameters (budgets, depth limit)
   !> @param[in]  seed_size     Leaf edge length (Bohr)
   !> @param[in]  refine_levels Refinement levels of a leaf without evidence
   !> @param[out] error         Invalid octree configuration
   !> @param[in]  debug         Per-anchor search report
   !> @param[in]  unit          Output unit of the report
   subroutine new_octree_seeder(seeder, param, seed_size, refine_levels, error, debug, unit)
      class(projection_seeder_type), allocatable, intent(out) :: seeder
      type(moist_cavity_drop_parameters_type), intent(in) :: param
      real(wp), intent(in) :: seed_size
      integer, intent(in) :: refine_levels
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in), optional :: debug
      integer, intent(in), optional :: unit

      allocate (octree_seeder_type :: seeder)
      select type (seeder)
      type is (octree_seeder_type)
         if (present(debug)) seeder%debug = debug
         if (present(unit)) seeder%unit = unit
         call seeder%init(error, seed_size=seed_size, max_boxes=param%octree_max_boxes, &
                          max_survivors=param%octree_max_survivors, max_depth=param%octree_max_depth, &
                          refine_levels=refine_levels)
      end select
   end subroutine new_octree_seeder

   !> Run the pipeline on one chunk of anchors
   !>
   !> @param[inout] self    Worker
   !> @param[in]    slots   Anchor slots, written into the rows (n)
   !> @param[in]    anchors Anchor positions (3, n)
   !> @param[in]    owners  Owning atoms (n)
   !> @param[in]    wlebs   Anchor quadrature weights (n)
   !> @param[in]    fs      Anchor switching values (n); below `wleb_cut` the anchor is skipped
   !> @param[out]   nrows   Rows appended per anchor (n)
   !> @param[out]   lost    Anchor kept no accepted foot (n)
   !> @param[out]   error   Build-fatal failure
   subroutine worker_process_chunk(self, slots, anchors, owners, wlebs, fs, nrows, lost, error)
      class(projection_worker_type), intent(inout) :: self
      integer, intent(in) :: slots(:)
      real(wp), intent(in) :: anchors(:, :)
      integer, intent(in) :: owners(:)
      real(wp), intent(in) :: wlebs(:)
      real(wp), intent(in) :: fs(:)
      integer, intent(out) :: nrows(:)
      logical, intent(out) :: lost(:)
      type(error_type), allocatable, intent(out) :: error

      integer :: k, r, j, jacc, icol, nnear, nnear_seeds
      logical :: accepted

      nrows = 0
      lost = .false.

      !> Pass 1: local Newton from every anchor above the switching cut
      call self%pass1%clear()
      do k = 1, size(slots)
         if (fs(k) < self%fcfg%wleb_cut) cycle
         call self%pass1%begin_run(k, anchors(:, k), owners(k))
         call self%pass1%append(anchors(:, k), 0.0_wp, icol)
      end do
      if (self%pass1%n == 0) return
      call newton_solve_batch(self%pass1, self%lsf, self%ncfg, self%nscratch, error)
      if (allocated(error)) return
      call polish_batch(self%pass1, self%lsf, self%pcfg, self%ncfg, self%nscratch, error)
      if (allocated(error)) return

      ! Seeding; a failed local solve waits for the nearest-foot search
      call self%pass2%clear()
      nnear = 0
      nnear_seeds = 0
      do r = 1, self%pass1%nrun
         k = self%pass1%run_anchor(r)
         j = self%pass1%run_start(r)
         accepted = run_has_accepted(self%pass1, r)
         if (.not. accepted .and. allocated(self%fallback)) then
            self%n_escalated = self%n_escalated + 1
            call queue_nearest(self%fallback, r, k)
         else if (.not. accepted .and. allocated(self%seeder)) then
            self%n_escalated = self%n_escalated + 1
            call queue_nearest(self%seeder, r, k)
         else if (allocated(self%seeder)) then
            jacc = 0
            if (self%pass1%status(j) == point_accepted) jacc = j
            call seed_anchor(self%seeder, jacc, r, k)
         else
            call finish_anchor(self%pass1, j, j, k)
         end if
         if (allocated(error)) return
      end do
      call solve_nearest()
      if (allocated(error)) return
      call flush_pass2()

   contains

      !> Request of the anchor of pass-1 run `r1`
      !>
      !> @param[in]  r1      Pass-1 run index
      !> @param[in]  k1      Anchor index within the chunk
      !> @param[out] request Request without a starting radius
      subroutine anchor_request(r1, k1, request)
         integer, intent(in) :: r1, k1
         type(seed_request_type), intent(out) :: request

         request%index = slots(k1)
         request%anchor = anchors(:, k1)
         request%owner_center = self%centers(:, owners(k1))
         call self%pass1%anchor_jet(r1, request%lsf0_anchor, request%grad_anchor)
         request%rho2_slack = self%fcfg%rho2_slack
         request%rho_fallback = fallback_search_radius
         request%rho_start = -1.0_wp
      end subroutine anchor_request

      !> Seed one pass-1 run and open its pass-2 run, flushing the batch first when full
      !>
      !> - the accepted local foot, if any, leads the run and bounds the search,
      !>   however far it lies
      !>
      !> @param[inout] seeder Seeder
      !> @param[in]    jacc   Pass-1 column of the accepted local foot, 0 for none
      !> @param[in]    r1     Pass-1 run index
      !> @param[in]    k1     Anchor index within the chunk
      subroutine seed_anchor(seeder, jacc, r1, k1)
         class(projection_seeder_type), intent(inout) :: seeder
         integer, intent(in) :: jacc, r1, k1
         type(seed_request_type) :: request
         integer :: n_new, i

         call anchor_request(r1, k1, request)
         if (jacc > 0) request%rho_start = sqrt(self%pass1%rho(jacc)**2 + request%rho2_slack)
         call seeder%seed(self%lsf, request, error)
         if (allocated(error)) return
         self%n_seeded = self%n_seeded + 1
         self%n_seeds = self%n_seeds + seeder%n_points

         n_new = seeder%n_points
         if (jacc > 0) n_new = n_new + 1
         if (self%pass2%would_overflow(n_new)) then
            call flush_pass2()
            if (allocated(error)) return
         end if
         call self%pass2%begin_run(k1, anchors(:, k1), owners(k1))
         if (jacc > 0) call self%pass2%append_copy(self%pass1, jacc, icol)
         ! Seeds start with a zero multiplier
         do i = 1, seeder%n_points
            call self%pass2%append(seeder%points(:, i), 0.0_wp, icol)
         end do
      end subroutine seed_anchor

      !> Store the nearest-first seeds of the anchor of pass-1 run `r1`
      !>
      !> - an anchor without a single seed is finished at once, as lost
      !>
      !> @param[inout] seeder Seeder of anchors without an accepted local foot
      !> @param[in]    r1     Pass-1 run index
      !> @param[in]    k1     Anchor index within the chunk
      subroutine queue_nearest(seeder, r1, k1)
         class(projection_seeder_type), intent(inout) :: seeder
         integer, intent(in) :: r1, k1
         type(seed_request_type) :: request
         integer :: i

         call anchor_request(r1, k1, request)
         request%nearest_first = .true.
         call seeder%seed(self%lsf, request, error)
         if (allocated(error)) return
         self%n_seeded = self%n_seeded + 1
         if (seeder%n_points == 0) then
            call finish_anchor(self%pass1, self%pass1%run_start(r1), self%pass1%run_start(r1), k1)
            return
         end if

         nnear = nnear + 1
         call reserve_nearest(nnear, nnear_seeds + seeder%n_points)
         self%near_run(nnear) = r1
         self%near_reach(nnear) = seeder%reach
         self%near_next(nnear) = nnear_seeds + 1
         do i = 1, seeder%n_points
            self%near_seeds(:, nnear_seeds + i) = seeder%points(:, i)
         end do
         nnear_seeds = nnear_seeds + seeder%n_points
         self%near_last(nnear) = nnear_seeds
      end subroutine queue_nearest

      !> Grow the nearest-foot queue to `n` anchors and `m` seeds, keeping its contents
      !>
      !> @param[in] n Anchors to hold
      !> @param[in] m Seeds to hold
      subroutine reserve_nearest(n, m)
         integer, intent(in) :: n, m
         integer, allocatable :: tmp(:)
         real(wp), allocatable :: tmp_reach(:), tmp_seeds(:, :)
         integer :: cap

         if (.not. allocated(self%near_run)) then
            allocate (self%near_run(64), self%near_last(64), self%near_next(64), self%near_reach(64))
            allocate (self%near_seeds(3, 1024))
         end if
         if (n > size(self%near_run)) then
            cap = max(n, 2*size(self%near_run))
            allocate (tmp_reach(cap))
            tmp_reach(1:n - 1) = self%near_reach(1:n - 1)
            call move_alloc(tmp_reach, self%near_reach)
            allocate (tmp(cap))
            tmp(1:n - 1) = self%near_run(1:n - 1)
            call move_alloc(tmp, self%near_run)
            allocate (tmp(cap))
            tmp(1:n - 1) = self%near_last(1:n - 1)
            call move_alloc(tmp, self%near_last)
            allocate (tmp(cap))
            tmp(1:n - 1) = self%near_next(1:n - 1)
            call move_alloc(tmp, self%near_next)
         end if
         if (m > size(self%near_seeds, 2)) then
            cap = max(m, 2*size(self%near_seeds, 2))
            allocate (tmp_seeds(3, cap))
            tmp_seeds(:, 1:nnear_seeds) = self%near_seeds(:, 1:nnear_seeds)
            call move_alloc(tmp_seeds, self%near_seeds)
         end if
      end subroutine reserve_nearest

      !> Nearest-foot search of every queued anchor, in rounds batched across anchors
      !>
      !> - round `i` solves the next `nearest_count 2^(i-1)` seeds of every
      !>   anchor still without an accepted foot
      !> - an anchor whose seeds run out without one is finished as lost
      subroutine solve_nearest()
         integer :: take, npend, nkeep, ip, e, r1, k1, rn, s, s_last, jbest, jn
         integer, allocatable :: pend(:)

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
            call self%passn%clear()
            do ip = 1, npend
               e = pend(ip)
               k1 = self%pass1%run_anchor(self%near_run(e))
               call self%passn%begin_run(e, anchors(:, k1), owners(k1))
               s_last = min(self%near_last(e), self%near_next(e) + take - 1)
               do s = self%near_next(e), s_last
                  call self%passn%append(self%near_seeds(:, s), 0.0_wp, icol)
               end do
               self%n_seeds = self%n_seeds + s_last - self%near_next(e) + 1
               self%near_next(e) = s_last + 1
            end do
            call newton_solve_batch(self%passn, self%lsf, self%ncfg, self%nscratch, error)
            if (allocated(error)) return
            call polish_batch(self%passn, self%lsf, self%pcfg, self%ncfg, self%nscratch, error)
            if (allocated(error)) return

            nkeep = 0
            do rn = 1, self%passn%nrun
               e = self%passn%run_anchor(rn)
               r1 = self%near_run(e)
               k1 = self%pass1%run_anchor(r1)
               jbest = 0
               do jn = self%passn%run_start(rn), self%passn%run_start(rn) + self%passn%run_len(rn) - 1
                  if (self%passn%status(jn) /= point_accepted) cycle
                  if (jbest == 0) then
                     jbest = jn
                  else if (self%passn%rho(jn) < self%passn%rho(jbest)) then
                     jbest = jn
                  end if
               end do
               if (jbest > 0) then
                  if (self%confirm_nearest) then
                     call confirm_nearest_foot(rn, e, jbest, k1)
                  else
                     call finish_anchor(self%passn, jbest, jbest, k1)
                  end if
                  if (allocated(error)) return
               else if (self%near_next(e) <= self%near_last(e)) then
                  nkeep = nkeep + 1
                  pend(nkeep) = e
               else
                  call finish_anchor(self%pass1, self%pass1%run_start(r1), self%pass1%run_start(r1), k1)
                  if (allocated(error)) return
               end if
            end do
            npend = nkeep
            take = 2*take
         end do
      end subroutine solve_nearest

      !> Open the pass-2 run of an anchor whose nearest-foot round accepted a foot
      !>
      !> - the round's accepted feet, then the unsolved seeds closer to the anchor
      !>   than `sqrt(rho^2 + slack)` plus the seeder's reach, `rho` of the
      !>   nearest accepted foot; the seeds are sorted, so the cut is a prefix
      !>
      !> @param[in] rn    Run of the anchor in the nearest-foot batch
      !> @param[in] e     Queue index of the anchor
      !> @param[in] jbest Nearest accepted column of the run
      !> @param[in] k1    Anchor index within the chunk
      subroutine confirm_nearest_foot(rn, e, jbest, k1)
         integer, intent(in) :: rn, e, jbest, k1
         real(wp) :: cut
         integer :: jn, j_first, j_last, s, s_last, n_new

         j_first = self%passn%run_start(rn)
         j_last = j_first + self%passn%run_len(rn) - 1
         cut = sqrt(self%passn%rho(jbest)**2 + self%fcfg%rho2_slack) + self%near_reach(e)
         s_last = self%near_next(e) - 1
         do s = self%near_next(e), self%near_last(e)
            if (norm2(self%near_seeds(:, s) - anchors(:, k1)) > cut) exit
            s_last = s
         end do

         n_new = s_last - self%near_next(e) + 1
         do jn = j_first, j_last
            if (self%passn%status(jn) == point_accepted) n_new = n_new + 1
         end do
         if (self%pass2%would_overflow(n_new)) then
            call flush_pass2()
            if (allocated(error)) return
         end if
         call self%pass2%begin_run(k1, anchors(:, k1), owners(k1))
         do jn = j_first, j_last
            if (self%passn%status(jn) == point_accepted) call self%pass2%append_copy(self%passn, jn, icol)
         end do
         do s = self%near_next(e), s_last
            call self%pass2%append(self%near_seeds(:, s), 0.0_wp, icol)
         end do
         self%n_seeds = self%n_seeds + s_last - self%near_next(e) + 1
      end subroutine confirm_nearest_foot

      !> Solve, polish and finalize every run of the pass-2 batch
      subroutine flush_pass2()
         integer :: r2, k2

         if (self%pass2%nrun == 0) return
         if (self%pass2%n > 0) then
            call newton_solve_batch(self%pass2, self%lsf, self%ncfg_seeded, self%nscratch, error)
            if (allocated(error)) return
            call polish_batch(self%pass2, self%lsf, self%pcfg, self%ncfg_seeded, self%nscratch, error)
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
         integer :: nbranch

         rec%index = slots(k2)
         rec%xyz = anchors(:, k2)
         rec%owner_center = self%centers(:, owners(k2))
         rec%wleb = wlebs(k2)
         rec%f = fs(k2)
         call finalize_anchor(self%fcfg, batch, first, last, rec, self%f_crit, self%f_foc, self%f_wleb, &
                              self%branch_weight, self%rows, nbranch, nrows(k2), self%fscratch, error)
         lost(k2) = nbranch == 0
      end subroutine finish_anchor

   end subroutine worker_process_chunk

   !> Project one anchor and return its rows (the single-anchor entry for tests)
   !>
   !> - the worker's row buffer is cleared first, so the outputs are exactly
   !>   this anchor's surviving branches in canonical order
   !>
   !> @param[inout] self      Worker
   !> @param[in]    anchor    Anchor position
   !> @param[in]    owner     Owning atom
   !> @param[out]   n_points  Rows produced (0: anchor lost or below the cut)
   !> @param[out]   error     Build-fatal failure
   !> @param[out]   feet      Projected points (3, n_points)
   !> @param[out]   normals   Unit normals (3, n_points)
   !> @param[out]   lambda    Multipliers
   !> @param[out]   rho       Anchor distances
   !> @param[out]   phi       Objective values
   !> @param[out]   weights   Branch weights
   !> @param[out]   n_seeds   Seeds the seeder proposed for this anchor
   !> @param[out]   lost      No accepted foot at all
   !> @param[in]    wleb      Anchor quadrature weight (default 1)
   !> @param[in]    f         Anchor switching value (default 1)
   subroutine worker_project_anchor(self, anchor, owner, n_points, error, feet, normals, lambda, rho, &
                                    phi, weights, n_seeds, lost, wleb, f)
      class(projection_worker_type), intent(inout) :: self
      real(wp), intent(in) :: anchor(3)
      integer, intent(in) :: owner
      integer, intent(out) :: n_points
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable, intent(out), optional :: feet(:, :), normals(:, :)
      real(wp), allocatable, intent(out), optional :: lambda(:), rho(:), phi(:), weights(:)
      integer, intent(out), optional :: n_seeds
      logical, intent(out), optional :: lost
      real(wp), intent(in), optional :: wleb
      real(wp), intent(in), optional :: f

      real(wp) :: wleb_k(1), f_k(1)
      integer :: nrows(1), seeds_before
      logical :: lost_k(1)

      wleb_k = 1.0_wp
      f_k = 1.0_wp
      if (present(wleb)) wleb_k = wleb
      if (present(f)) f_k = f

      call self%rows%clear()
      seeds_before = self%n_seeds
      call self%process_chunk([1], reshape(anchor, [3, 1]), [owner], wleb_k, f_k, nrows, lost_k, error)
      n_points = 0
      if (allocated(error)) return
      n_points = self%rows%n
      if (present(n_seeds)) n_seeds = self%n_seeds - seeds_before
      if (present(lost)) lost = lost_k(1)
      if (present(feet)) feet = self%rows%xyz(:, 1:n_points)
      if (present(normals)) normals = self%rows%normal(:, 1:n_points)
      if (present(lambda)) lambda = self%rows%lambda(1:n_points)
      if (present(rho)) rho = self%rows%rho(1:n_points)
      if (present(phi)) phi = self%rows%phi(1:n_points)
      if (present(weights)) weights = self%rows%wbranch(1:n_points)
   end subroutine worker_project_anchor

   !> Whether any column of a run was accepted
   !>
   !> @param[in] batch Batch
   !> @param[in] irun  Run index
   pure function run_has_accepted(batch, irun) result(found)
      !> Batch
      type(projection_batch_type), intent(in) :: batch
      !> Run index
      integer, intent(in) :: irun
      !> Whether an accepted column exists
      logical :: found

      integer :: j

      found = .false.
      do j = batch%run_start(irun), batch%run_start(irun) + batch%run_len(irun) - 1
         if (batch%status(j) == point_accepted) then
            found = .true.
            return
         end if
      end do
   end function run_has_accepted

end module moist_cavity_drop_projection_worker
