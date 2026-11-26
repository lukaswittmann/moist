!> Optimality verification and saddle escape for converged DROP projections
!>
!> Per-point routines ported from the single-anchor projector:
!>   - [[check_constrained_optimality]] -- KKT residuals plus the reduced
!>     Hessian on the tangent plane, third order when that is flat
!>   - [[riemannian_newton_escape]] -- walk off a saddle along the surface
!>   - [[retract_to_surface]] -- Newton along the normal back onto `S = 0`
!>
!> [[polish_batch]] applies them to every converged column of a batch, reusing
!> the jet Newton left there, and re-polishes escaped saddles, constrained maxima
!> included, with batched Newton passes
!>
!> Run context read only, shared by the worker threads: per-anchor traces through
!> `ctx%message` at verbosity 6 (`debug` raises to it), per-iteration detail
!> through `ctx%debug_message`; no timers or warnings
!>
!> Every level set evaluation goes through the LSF's `evaluate_batch`, the third
!> derivative of the degenerate test included; a point beyond the level set
!> reach is a recoverable failure of the routine that met it, never build-fatal
module moist_cavity_drop_projection_polish
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_math_packed, only: packed_full_contract3
   use moist_utils_format, only: to_string
   use moist_utils_mem, only: grow_array
   use moist_cavity_drop_projection_batch, only: projection_batch_type, &
                                                point_pending, point_converged, point_accepted, point_captured, &
                                                point_failed
   use moist_cavity_drop_projection_newton, only: newton_scratch_type, newton_solve_batch, captor_list_type, &
                                                  kkt_residuals, reduced_hessian
   use moist_cavity_drop_projection_parameters, only: drop_solve_parameters_type, drop_polish_parameters_type
   implicit none(type, external)
   private

   public :: polish_batch, saddle_cache_type
   public :: check_constrained_optimality, riemannian_newton_escape, retract_to_surface

   !> Third directional derivative of the Lagrangian above which a flat point is an inflection
   real(wp), parameter :: tol_D3 = 1.0e-8_wp
   !> Largest objective drop an escape step must achieve; less is ok for a step with a smaller
   !> curvature-based gain at the threshold
   real(wp), parameter :: escape_drop_max = 1.0e-10_wp
   !> Fixed step along negative curvature when leaving a saddle
   real(wp), parameter :: escape_step = 0.5_wp
   !> Line-search trials per escape iteration
   integer, parameter :: escape_max_ls_iter = 20
   !> Reduction of a rejected escape line-search step
   real(wp), parameter :: escape_alpha_reduce = 0.75_wp
   !> Armijo coefficient of the escape line search
   real(wp), parameter :: escape_armijo_c = 1.0e-4_wp
   !> Squared distance (Bohr^2) within which two saddles of one anchor match for deduplication
   real(wp), parameter :: saddle_dedup_tol2 = 1.0e-12_wp
   !> Relative eigenvalue gap below which an index-2 saddle is isotropic and never deduplicated
   real(wp), parameter :: saddle_iso_rel = 1.0e-3_wp
   !> Relative objective gap below which a saddle's two escape signs are near-symmetric
   real(wp), parameter :: saddle_sym_rel = 1.0e-3_wp

   !> State of a saddle escape direction: untaken, walk not yet settled, confirmed, or dropped entry
   integer, parameter :: lead_free = 0, lead_pending = 1, lead_confirmed = 2, lead_dropped = 3
   !> A visit of a cached saddle: captured, waiting for its leader, or walking the opposite sign
   integer, parameter :: visit_capture = 1, visit_wait = 2, visit_opposite = 3

   !> Per-anchor cache of escaped saddles for symmetry-safe duplicate skipping
   !>
   !> - keyed by an anchor key (`1..nkey`): the run index within one batch, or the
   !>   run's `run_key`, persistent across the rounds of the extended search
   !> - one entry per non-isotropic saddle: its position, the unit 3D escape
   !>   direction of its first walk, whether the two signs were near-symmetric,
   !>   and per direction (first, opposite) its state: free, led by a walk not yet
   !>   settled, or confirmed
   !> - confirmed: the walk ended `point_accepted` or captured by an accepted point,
   !>   so the point it leads to is in the run; a walk that fails drops the entry
   !> - [[polish_batch]] settles every walk before it returns, so the entries a
   !>   persistent cache carries across calls hold free and confirmed directions only
   type :: saddle_cache_type
      !> Active anchor keys (`1:nkey`)
      integer :: nkey = 0
      !> Linked-list head saddle index per anchor key (0 for none)
      integer, allocatable :: head(:)
      !> Stored saddles across active keys, dropped ones included
      integer :: n = 0
      !> Next saddle index in the same anchor key's list (0 for end)
      integer, allocatable :: next(:)
      !> Anchor key of the saddle
      integer, allocatable :: key(:)
      !> Saddle positions (3, :) and unit 3D escape direction of the first walk (3, :)
      real(wp), allocatable :: x(:, :), dir(:, :)
      !> Whether the saddle's initial escape signs had near-equal objectives; cleared
      !> once the opposite sign showed no descent
      logical, allocatable :: near_sym(:)
      !> State of the first and the opposite direction (2, :): `lead_free`,
      !> `lead_pending`, `lead_confirmed` or `lead_dropped`
      integer, allocatable :: lead(:, :)
   contains
      !> Clear all stored saddles while preserving allocated storage
      procedure :: reset => saddle_cache_reset
      !> Release all storage
      procedure :: destroy => saddle_cache_destroy
   end type saddle_cache_type

contains

   !* ================================================================================= *!
   !*                                    Batch polish                                   *!
   !* ================================================================================= *!

   !> Verify every converged column in `first:last`, escaping saddles
   !>
   !> - minimum: `point_accepted`; no KKT point: `point_failed`
   !> - saddle (a constrained maximum included): escaped, re-solved by a batched
   !>   Newton pass and re-checked; a point that is still a saddle afterwards is
   !>   `point_failed`
   !> - a classification the third-order test cannot settle keeps the point,
   !>   as the single-anchor projector did
   !>
   !> With `capture_duplicates` (of `ncfg`), per anchor:
   !>   - the escape walk stops as `point_captured` once it lands on an accepted
   !>     point of its run
   !>   - a repeated visit of a non-isotropic saddle is `point_captured` without a
   !>     walk once every direction it could take is confirmed (one direction, both
   !>     signs of a near-symmetric saddle, whose second visit walks the opposite
   !>     sign); a visit whose leader is not settled waits for the re-check of its
   !>     round, a failed leader drops the entry and a waiting visit walks itself
   !>   - an opposite sign that finds no descent at the saddle (an inflection, a
   !>     saddle its third order dominates) is no failed walk: the entry stays,
   !>     one-sided, and the visit is a repeat of the first direction
   !>   - rounds repeat until no visit waits; each round settles every walk of
   !>     the previous one, so its first visit of a saddle never waits, and the
   !>     outcome depends on the column order of the batch only
   !>
   !> @param[in]    ctx          Run context, read only; traces at verbosity 6, detail on the debug stream
   !> @param[inout] batch        Points with Newton's final jets
   !> @param[inout] lsf          Thread-local level set
   !> @param[in]    cfg          Polish settings
   !> @param[in]    ncfg         Newton settings for the re-polish
   !> @param[in]    phi_alpha    Anchor objective weight
   !> @param[inout] nscratch     Newton scratch
   !> @param[out]   error        LSF evaluation failure (build-fatal)
   !> @param[in]    first        First column (default 1)
   !> @param[in]    last         Last column (default `batch%n`)
   !> @param[inout] saddle_cache Persistent per-anchor saddle cache (keyed by `batch%run_key`)
   subroutine polish_batch(ctx, batch, lsf, cfg, ncfg, phi_alpha, nscratch, error, first, last, saddle_cache)
      !> Run context, read only; traces at verbosity 6, detail on the debug stream
      type(moist_context_type), intent(in) :: ctx
      !> Points with Newton's final jets
      type(projection_batch_type), intent(inout) :: batch
      !> Thread-local level set
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Polish settings
      type(drop_polish_parameters_type), intent(in) :: cfg
      !> Newton settings for the re-polish
      type(drop_solve_parameters_type), intent(in) :: ncfg
      !> Anchor objective weight
      real(wp), intent(in) :: phi_alpha
      !> Newton scratch
      type(newton_scratch_type), intent(inout) :: nscratch
      !> LSF evaluation failure (build-fatal)
      type(error_type), allocatable, intent(out) :: error
      !> First column (default 1)
      integer, intent(in), optional :: first
      !> Last column (default `batch%n`)
      integer, intent(in), optional :: last
      !> Persistent per-anchor saddle cache (keyed by `batch%run_key`)
      type(saddle_cache_type), intent(inout), optional, target :: saddle_cache

      type(error_type), allocatable :: escape_error
      type(saddle_cache_type), target :: local_cache
      type(saddle_cache_type), pointer :: scache
      !> Saddle columns in column order and whether each is isotropic
      integer, allocatable :: sad_col(:)
      logical, allocatable :: sad_iso(:)
      !> Cache entry and direction (1 first, 2 opposite) each saddle leads, 0 for none
      integer, allocatable :: sad_entry(:), sad_sign(:)
      !> Saddles of the round, saddles waiting for the next one, saddles walked in the round
      integer, allocatable :: todo(:), waiting(:), walked(:)
      !> Direction the walk must not take; run whose captors stop the walk (unallocated: absent)
      real(wp), allocatable :: avoid(:)
      integer, allocatable :: cap_run
      real(wp) :: r_out(3), lambda_out, mu_min, mu_max, esc_dir(3)
      integer :: j0, j1, j, k, i, irun, skey, m, n_sad, n_todo, n_wait, n_walk
      logical :: saddle, kkt, do_cap, is_sym, esc_cap

      j0 = 1
      j1 = batch%n
      if (present(first)) j0 = first
      if (present(last)) j1 = last
      if (j1 < j0) return
      do_cap = ncfg%capture_duplicates
      allocate (sad_col(j1 - j0 + 1), sad_iso(j1 - j0 + 1))

      if (do_cap) then
         call nscratch%captors%init(batch, batch%run_start(batch%run_of(j0)), j1)
         if (present(saddle_cache)) then
            scache => saddle_cache
         else
            scache => local_cache
            call scache%reset(batch%nrun)
         end if
      else
         nullify (scache)
      end if

      !> Pass 1: classify every converged column first so all non-saddle minima in
      !> the batch are registered as captors before any saddle escape walk begins
      n_sad = 0
      do j = j0, j1
         if (batch%status(j) /= point_converged) cycle
         call classify(j, saddle, kkt, mu_min, mu_max)
         if (allocated(error)) return
         if (.not. kkt) then
            batch%status(j) = point_failed
         else if (.not. saddle) then
            batch%status(j) = point_accepted
            if (do_cap) call nscratch%captors%add(batch, j)
         else
            n_sad = n_sad + 1
            sad_col(n_sad) = j
            sad_iso(n_sad) = (mu_max - mu_min) <= saddle_iso_rel*max(1.0_wp, abs(mu_min))
         end if
      end do
      if (n_sad == 0) return

      !> Pass 2: escape the saddles in column order, in rounds
      allocate (sad_entry(n_sad), sad_sign(n_sad), todo(n_sad), waiting(n_sad), walked(n_sad))
      do k = 1, n_sad
         todo(k) = k
      end do
      n_todo = n_sad
      do while (n_todo > 0)
         n_wait = 0
         n_walk = 0
         do i = 1, n_todo
            k = todo(i)
            j = sad_col(k)
            irun = batch%run_of(j)
            skey = irun
            if (present(saddle_cache)) skey = batch%run_key(irun)
            sad_entry(k) = 0
            sad_sign(k) = 1
            if (allocated(avoid)) deallocate (avoid)

            !> Symmetry-safe per-anchor saddle deduplication
            if (do_cap .and. .not. sad_iso(k)) then
               m = saddle_cache_find(scache, skey, batch%x(:, j))
               if (m > 0) then
                  select case (saddle_cache_visit(scache, m))
                  case (visit_capture)
                     batch%status(j) = point_captured
                     cycle
                  case (visit_wait)
                     n_wait = n_wait + 1
                     waiting(n_wait) = k
                     cycle
                  case default
                     sad_entry(k) = m
                     sad_sign(k) = 2
                     scache%lead(2, m) = lead_pending
                     avoid = scache%dir(:, m)
                  end select
               end if
            end if

            if (ctx%verbosity(6)) then
               call ctx%message("[optimality] anchor "//to_string(batch%run_id(irun))// &
                                ": Saddle point detected, attempting escape ...", level=6)
            end if
            if (do_cap) cap_run = irun
            call riemannian_newton_escape(ctx, cfg, phi_alpha, ncfg%max_step, lsf, batch%x(:, j), &
                                          batch%run_xyz(:, irun), batch%lambda(j), r_out, lambda_out, &
                                          escape_error, lsf_error=error, jet_S=batch%S(j), &
                                          jet_grad=batch%gradS(:, j), jet_hess=batch%hessS(:, j), &
                                          avoid_dir=avoid, escape_dir=esc_dir, is_near_sym=is_sym, &
                                          captor_batch=batch, captors=nscratch%captors, capture_run=cap_run, &
                                          captured=esc_cap)
            if (allocated(error)) return
            ! A walk that ran out of iterations still hands back its best point;
            ! Newton and the re-check below decide whether it is usable
            if (allocated(escape_error)) deallocate (escape_error)

            ! An opposite sign without descent leads nowhere: the saddle is one-sided
            ! from here on and the visit a plain repeat of its first direction
            if (sad_sign(k) == 2 .and. all(esc_dir == 0.0_wp)) then
               m = sad_entry(k)
               scache%near_sym(m) = .false.
               scache%lead(2, m) = lead_free
               if (saddle_cache_visit(scache, m) == visit_capture) then
                  batch%status(j) = point_captured
               else
                  n_wait = n_wait + 1
                  waiting(n_wait) = k
               end if
               cycle
            end if

            ! A first visit opens the saddle's entry and leads its first direction
            if (do_cap .and. .not. sad_iso(k) .and. sad_entry(k) == 0 .and. &
                dot_product(esc_dir, esc_dir) > 0.0_wp) then
               call saddle_cache_add(scache, skey, batch%x(:, j), esc_dir, is_sym, sad_entry(k))
            end if

            if (esc_cap) then
               ! Stopped on an accepted point: the direction is confirmed
               batch%status(j) = point_captured
               if (sad_entry(k) > 0) call saddle_cache_settle(scache, sad_entry(k), sad_sign(k), .true.)
               cycle
            end if

            batch%x(:, j) = r_out
            batch%lambda(j) = lambda_out
            batch%status(j) = point_pending
            n_walk = n_walk + 1
            walked(n_walk) = k
         end do

         if (n_walk > 0) then
            call newton_solve_batch(batch, lsf, ncfg, cfg, phi_alpha, nscratch, error, first=j0, last=j1)
            if (allocated(error)) return
            do i = 1, n_walk
               k = walked(i)
               j = sad_col(k)
               if (batch%status(j) == point_converged) then
                  call classify(j, saddle, kkt, mu_min, mu_max)
                  if (allocated(error)) return
                  if (.not. kkt) then
                     batch%status(j) = point_failed
                  else if (saddle) then
                     if (ctx%verbosity(6)) then
                        call ctx%message("[optimality] anchor "//to_string(batch%run_id(batch%run_of(j)))// &
                                         ": Failed to escape saddle point!", level=6)
                     end if
                     batch%status(j) = point_failed
                  else
                     batch%status(j) = point_accepted
                  end if
               end if
               ! Confirmed: accepted, or captured by Newton on a captor the check accepts
               if (sad_entry(k) > 0) then
                  call saddle_cache_settle(scache, sad_entry(k), sad_sign(k), &
                                           batch%status(j) == point_accepted .or. batch%status(j) == point_captured)
               end if
            end do
         end if

         ! Waiting visits form the next round; the captors of its walks are the
         ! points accepted so far
         todo(1:n_wait) = waiting(1:n_wait)
         n_todo = n_wait
         if (n_todo > 0) call nscratch%captors%init(batch, batch%run_start(batch%run_of(j0)), j1)
      end do

   contains

      !> Classify converged column `j`; an unreliable classification is traced
      !>
      !> @param[in]  j      Column
      !> @param[out] saddle Saddle or inflection
      !> @param[out] kkt    Feasible and stationary
      !> @param[out] mu_min Smaller reduced-Hessian eigenvalue
      !> @param[out] mu_max Larger reduced-Hessian eigenvalue
      subroutine classify(j, saddle, kkt, mu_min, mu_max)
         integer, intent(in) :: j
         logical, intent(out) :: saddle
         logical, intent(out) :: kkt
         real(wp), intent(out) :: mu_min
         real(wp), intent(out) :: mu_max

         type(error_type), allocatable :: class_error
         integer :: irun

         irun = batch%run_of(j)
         call check_constrained_optimality(ctx, cfg, phi_alpha, lsf, batch%x(:, j), batch%run_xyz(:, irun), &
                                           batch%lambda(j), class_error, is_saddle=saddle, is_kkt=kkt, &
                                           anchor_id=batch%run_id(irun), jet_S=batch%S(j), &
                                           jet_grad=batch%gradS(:, j), jet_hess=batch%hessS(:, j), &
                                           lsf_error=error, mu_min_out=mu_min, mu_max_out=mu_max)
         if (.not. allocated(class_error)) return
         if (ctx%verbosity(6)) then
            call ctx%message("[optimality] anchor "//to_string(batch%run_id(irun))//": "// &
                             trim(class_error%message), level=6)
         end if
      end subroutine classify

   end subroutine polish_batch

   !* ================================================================================= *!
   !*                            Check constrained optimality                           *!
   !* ================================================================================= *!

   !> Whether a converged KKT point is a constrained local minimum or a saddle
   !>
   !> - reduced Hessian eigenvalues on the tangent plane, not the bordered KKT
   !>   matrix
   !> - a constrained maximum is a saddle here: both eigenvalues are negative
   !> - a flat second-order test falls through to `D^3 L[d,d,d]` along the flat
   !>   direction; nonzero means an inflection, treated as an escapable saddle
   !> - flat to third order raises `error`, classification unreliable
   !> - `r_star` beyond the level set reach: `is_kkt = .false.` and `error`
   !> - a point, multiplier or jet that is not finite: `is_kkt = .false.` and `error`
   !>
   !> @param[in]    ctx          Run context, read only
   !> @param[in]    cfg          Polish settings
   !> @param[in]    phi_alpha    Anchor objective weight
   !> @param[inout] lsf          Thread-local level set
   !> @param[in]    r_star       Candidate point
   !> @param[in]    anchor       Anchor
   !> @param[in]    lambda       Multiplier at `r_star`
   !> @param[out]   error        Unreliable classification, or saddle without `is_saddle`
   !> @param[out]   is_saddle    Saddle or inflection
   !> @param[out]   is_kkt       Feasible and stationary
   !> @param[out]   is_converged Strict local minimum
   !> @param[in]    anchor_id    Global anchor ID for traces
   !> @param[in]    jet_S        Level set value at `r_star`, skips the evaluation
   !> @param[in]    jet_grad     Level set gradient at `r_star`
   !> @param[in]    jet_hess     Level set Hessian at `r_star`, packed (6)
   !> @param[out]   lsf_error    LSF evaluation failure; folded into `error` when absent
   !> @param[out]   mu_min_out   Smaller reduced-Hessian eigenvalue
   !> @param[out]   mu_max_out   Larger reduced-Hessian eigenvalue
   subroutine check_constrained_optimality(ctx, cfg, phi_alpha, lsf, r_star, anchor, lambda, error, &
                                           is_saddle, is_kkt, is_converged, anchor_id, &
                                           jet_S, jet_grad, jet_hess, lsf_error, &
                                           mu_min_out, mu_max_out)
      type(moist_context_type), intent(in) :: ctx
      type(drop_polish_parameters_type), intent(in) :: cfg
      real(wp), intent(in) :: phi_alpha
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: r_star(3)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: lambda
      type(error_type), allocatable, intent(out) :: error
      logical, intent(out), optional :: is_saddle
      logical, intent(out), optional :: is_kkt
      logical, intent(out), optional :: is_converged
      integer, intent(in), optional :: anchor_id
      real(wp), intent(in), optional :: jet_S
      real(wp), intent(in), optional :: jet_grad(3)
      real(wp), intent(in), optional :: jet_hess(:)
      type(error_type), allocatable, intent(out), optional :: lsf_error
      real(wp), intent(out), optional :: mu_min_out
      real(wp), intent(out), optional :: mu_max_out

      type(error_type), allocatable :: eval_error
      real(wp) :: S_val, grad_S(3), hess_S(6)
      real(wp) :: norm_grad_S, res_feas, res_stat
      real(wp) :: T_mat(3, 2), t1(3), t2(3)
      real(wp) :: mu_min, mu_max, v_min(2), v_max(2), eig_tol
      real(wp) :: d_degen(3), D3_ddd
      integer :: id
      logical :: unreached

      if (present(is_saddle)) is_saddle = .false.
      if (present(is_kkt)) is_kkt = .true.
      if (present(is_converged)) is_converged = .false.
      if (present(mu_min_out)) mu_min_out = 0.0_wp
      if (present(mu_max_out)) mu_max_out = 0.0_wp
      id = 0
      if (present(anchor_id)) id = anchor_id

      if (present(jet_S) .and. present(jet_grad) .and. present(jet_hess)) then
         S_val = jet_S
         grad_S = jet_grad
         hess_S = jet_hess
      else
         call point_jet(lsf, r_star, 2, S_val, unreached, eval_error, grad_S, hess_S)
         if (allocated(eval_error)) then
            call hand_over(eval_error, error, lsf_error)
            if (present(is_kkt)) is_kkt = .false.
            return
         end if
         if (unreached) then
            call fatal_error(error, "KKT check: candidate point beyond the level set reach")
            if (present(is_kkt)) is_kkt = .false.
            return
         end if
      end if

      call kkt_residuals(phi_alpha, r_star, anchor, lambda, S_val, grad_S, res_feas, res_stat)
      ! NaN fails every comparison below and would end as a kept point; finite
      ! residuals mean a finite point, multiplier, value and gradient
      if (.not. (ieee_is_finite(res_feas) .and. ieee_is_finite(res_stat) .and. all(ieee_is_finite(hess_S)))) then
         call fatal_error(error, "KKT check: candidate point or its level set jet not finite")
         if (present(is_kkt)) is_kkt = .false.
         return
      end if
      norm_grad_S = norm2(grad_S)
      if (norm_grad_S < cfg%kkt_tol_grad) then
         if (ctx%verbosity(6)) then
            call ctx%message("[KKT check] anchor "//to_string(id)//": Degenerate LSF gradient ||grad_S|| = "// &
                             to_string(norm_grad_S, "(es12.4)"), level=6)
         end if
         ! Cannot classify, but an infeasible point is no KKT point
         if (present(is_kkt)) is_kkt = res_feas <= cfg%kkt_tol_feas
         return
      end if

      if (ctx%debug) then
         ! One message, so no other thread's line lands inside it
         call ctx%debug_message(" KKT verification, anchor "//to_string(id)//":"//new_line("a")// &
                                "   Feasibility   |S(r*)| = "//to_string(res_feas, "(es12.4)")// &
                                "  (tol: "//to_string(cfg%kkt_tol_feas, "(es12.4)")//")"//new_line("a")// &
                                "   Stationarity  ||dL||  = "//to_string(res_stat, "(es12.4)")// &
                                "  (tol: "//to_string(cfg%kkt_tol_stat, "(es12.4)")//")")
      end if

      if (res_feas > cfg%kkt_tol_feas .or. res_stat > cfg%kkt_tol_stat) then
         if (ctx%verbosity(6)) then
            call ctx%message("[KKT check] anchor "//to_string(id)//": Not at KKT point: |S| = "// &
                             to_string(res_feas, "(es12.4)")//", ||dL|| = "//to_string(res_stat, "(es12.4)"), &
                             level=6)
         end if
         if (present(is_kkt)) is_kkt = .false.
         return
      end if

      !> Reduced Hessian R = T^T (Hess phi - lambda Hess S) T on the tangent plane
      call reduced_hessian(phi_alpha, lambda, grad_S, hess_S, cfg%saddle_eig_tol, t1, t2, &
                           mu_min, mu_max, v_min, v_max, eig_tol)
      T_mat(:, 1) = t1
      T_mat(:, 2) = t2
      if (present(mu_min_out)) mu_min_out = mu_min
      if (present(mu_max_out)) mu_max_out = mu_max

      if (ctx%debug) then
         call ctx%debug_message("   Reduced Hessian eigenvalues, anchor "//to_string(id)//": "// &
                                to_string(mu_min, "(es14.6)")//" "//to_string(mu_max, "(es14.6)"))
      end if

      if (mu_min < -eig_tol) then
         if (present(is_saddle)) then
            is_saddle = .true.
         else
            call fatal_error(error, "Saddle point detected but strict minimum required")
         end if
         return
      end if

      if (mu_min > eig_tol) then
         if (present(is_converged)) is_converged = .true.
         return
      end if

      !> Flat second order: D^3 L[d,d,d] = -lambda S_ijk d_i d_j d_k along d = T v_min
      d_degen = matmul(T_mat, v_min)
      call lagrangian_third_derivative(lsf, r_star, lambda, d_degen, D3_ddd, unreached, eval_error)
      if (allocated(eval_error)) then
         call hand_over(eval_error, error, lsf_error)
         if (present(is_kkt)) is_kkt = .false.
         return
      end if
      if (unreached) then
         call fatal_error(error, "KKT check: candidate point beyond the level set reach")
         if (present(is_kkt)) is_kkt = .false.
         return
      end if

      if (ctx%debug) then
         call ctx%debug_message("   Reduced eigs, anchor "//to_string(id)//" = "//to_string(mu_min, "(es16.8)")// &
                                " "//to_string(mu_max, "(es16.8)")//", D3 L[d,d,d] = "//to_string(D3_ddd, "(es16.8)"))
      end if

      if (abs(D3_ddd) > tol_D3) then
         if (present(is_saddle)) then
            is_saddle = .true.
         else
            call fatal_error(error, "Inflection point detected but strict minimum required")
         end if
      else
         call fatal_error(error, "Degenerate to third order; classification unreliable")
         ! Kept, unless the third derivative is NaN
         if (present(is_kkt)) is_kkt = abs(D3_ddd) <= tol_D3
      end if
   end subroutine check_constrained_optimality

   !* ================================================================================= *!
   !*                    Riemannian Newton escape from saddle points                    *!
   !* ================================================================================= *!

   !> Walk off a saddle along the surface, then descend to a minimum
   !>
   !> - at a saddle (a constrained maximum included), or at a flat point whose
   !>   third derivative makes it an inflection (the classification of
   !>   [[check_constrained_optimality]]): fixed step along the negative-curvature
   !>   or flat tangent direction; the first-order change there is roundoff at a
   !>   KKT point, so both signs are tried and the lower objective kept
   !> - a flat point with zero third derivative ends the walk, as the check keeps it
   !> - elsewhere: Riemannian Newton (or steepest descent when nearly singular),
   !>   the step at most `max_step` long
   !> - stops in the basin of a minimum, at a tangential gradient of `1e-6`;
   !>   [[polish_batch]] converges the point with the batched Newton like any other
   !>   column
   !> - every trial is retracted onto `S = 0` before the acceptance test; an escape
   !>   step must lower `phi`, a Newton or descent step must pass the Armijo test
   !> - `r_out`/`lambda_out` always hold the last accepted point, also on error
   !> - an iterate beyond the level set reach ends the walk with `error`
   !> - with `capture_run`: stops as `captured` at a converged point, or before an
   !>   uncapped Newton step on a convex patch, that lands on a captor of the run
   !> - `iter = 0..escape_max_iter`: a first step, then up to `escape_max_iter`
   !>   more
   !>
   !> @param[in]    ctx          Run context, read only; iterations on the debug stream
   !> @param[in]    cfg          Polish settings
   !> @param[in]    phi_alpha    Anchor objective weight
   !> @param[in]    max_step     Longest Newton step (Bohr), the Newton `max_step`
   !> @param[inout] lsf          Thread-local level set
   !> @param[in]    r_init       Start point
   !> @param[in]    anchor       Anchor
   !> @param[in]    lambda_init  Multiplier at `r_init`
   !> @param[out]   r_out        Final point
   !> @param[out]   lambda_out   Multiplier at `r_out`
   !> @param[out]   error        Degenerate gradient, unreached iterate or iteration budget exhausted
   !> @param[out]   lsf_error    LSF evaluation failure; folded into `error` when absent
   !> @param[in]    jet_S        Level set value at `r_init`, skips the `iter = 0` LSF call
   !> @param[in]    jet_grad     Level set gradient at `r_init`
   !> @param[in]    jet_hess     Level set Hessian at `r_init`, packed (6)
   !> @param[in]    avoid_dir    3D unit escape direction to avoid at `iter = 0` (forces the opposite sign)
   !> @param[out]   escape_dir   Unit 3D escape direction chosen at `iter = 0` (zero if none)
   !> @param[out]   is_near_sym  Initial escape signs had near-equal retracted objectives
   !> @param[in]    captor_batch Batch holding the captor columns of `capture_run`
   !> @param[in]    captors      Captor lists of `captor_batch`
   !> @param[in]    capture_run  Run whose captors stop the walk; capture is on only when present
   !> @param[out]   captured     The walk stopped on a captor
   subroutine riemannian_newton_escape(ctx, cfg, phi_alpha, max_step, lsf, r_init, anchor, lambda_init, r_out, &
                                       lambda_out, error, lsf_error, jet_S, jet_grad, jet_hess, avoid_dir, &
                                       escape_dir, is_near_sym, captor_batch, captors, capture_run, captured)
      type(moist_context_type), intent(in) :: ctx
      type(drop_polish_parameters_type), intent(in) :: cfg
      real(wp), intent(in) :: phi_alpha
      real(wp), intent(in) :: max_step
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: r_init(3)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: lambda_init
      real(wp), intent(out) :: r_out(3)
      real(wp), intent(out) :: lambda_out
      type(error_type), allocatable, intent(out) :: error
      type(error_type), allocatable, intent(out), optional :: lsf_error
      real(wp), intent(in), optional :: jet_S
      real(wp), intent(in), optional :: jet_grad(3)
      real(wp), intent(in), optional :: jet_hess(:)
      real(wp), intent(in), optional :: avoid_dir(3)
      real(wp), intent(out), optional :: escape_dir(3)
      logical, intent(out), optional :: is_near_sym
      type(projection_batch_type), intent(in), optional :: captor_batch
      type(captor_list_type), intent(in), optional :: captors
      integer, intent(in), optional :: capture_run
      logical, intent(out), optional :: captured

      !> Hand-over gradient: below about sqrt(eps) the objective-based line search
      !> cannot resolve progress, and a point handed over below Newton's residual
      !> tolerance would be kept unrefined, off by `||F|| / mu` on a soft point
      real(wp), parameter :: tol_grad_tan = 1.0e-6_wp

      type(error_type), allocatable :: eval_error, retract_error
      real(wp) :: r_curr(3), r_trial(3), r_retracted(3), lambda_curr, lambda_trial
      real(wp) :: S_val, grad_S(3), hess_S(6), grad_phi(3)
      real(wp) :: phi_init, phi_curr, phi_trial, norm_grad_S
      real(wp) :: T_mat(3, 2), t1(3), t2(3), g_tan(2)
      real(wp) :: delta_r(3), alpha, directional_deriv, norm_g_tan
      real(wp) :: mu_min, mu_max, v_min(2), v_max(2), s_tan(2), eig_tol, D3_ddd
      real(wp) :: r_cand(3), lambda_cand, phi_cand, phi_pos, phi_neg, drop_min, dlen
      integer :: iter, ls_iter, isign, nsign, chosen_sign
      logical :: converged, ls_success, at_saddle, escaping, have_trial, unreached
      logical :: have_jet0, have_pos, have_neg, check_cap

      have_jet0 = present(jet_S) .and. present(jet_grad) .and. present(jet_hess)
      check_cap = present(captor_batch) .and. present(captors) .and. present(capture_run) .and. present(captured)
      if (present(escape_dir)) escape_dir = 0.0_wp
      if (present(is_near_sym)) is_near_sym = .false.
      if (present(captured)) captured = .false.
      r_curr = r_init
      lambda_curr = lambda_init
      r_out = r_curr
      lambda_out = lambda_curr
      converged = .false.
      phi_curr = 0.5_wp*phi_alpha*sum((r_curr - anchor)**2)
      phi_init = phi_curr

      do iter = 0, cfg%escape_max_iter
         if (iter == 0 .and. have_jet0) then
            S_val = jet_S
            grad_S = jet_grad
            hess_S = jet_hess
            unreached = .false.
         else
            call point_jet(lsf, r_curr, 2, S_val, unreached, eval_error, grad_S, hess_S)
            if (allocated(eval_error)) then
               call hand_over(eval_error, error, lsf_error)
               return
            end if
            if (unreached) then
               call fatal_error(error, "Riemannian escape: iterate beyond the level set reach")
               return
            end if
         end if
         grad_phi = phi_alpha*(r_curr - anchor)

         norm_grad_S = norm2(grad_S)
         if (norm_grad_S < cfg%kkt_tol_grad) then
            call fatal_error(error, "Riemannian escape: degenerate LSF gradient")
            return
         end if
         ! Same reduced Hessian and threshold as the optimality check, so both agree on what a saddle is
         call reduced_hessian(phi_alpha, lambda_curr, grad_S, hess_S, cfg%saddle_eig_tol, t1, t2, &
                              mu_min, mu_max, v_min, v_max, eig_tol)
         T_mat(:, 1) = t1
         T_mat(:, 2) = t2

         g_tan(1) = dot_product(t1, grad_phi)
         g_tan(2) = dot_product(t2, grad_phi)
         norm_g_tan = norm2(g_tan)
         at_saddle = mu_min < -eig_tol

         escaping = at_saddle
         if (norm_g_tan < tol_grad_tan .and. .not. at_saddle) then
            if (mu_min > eig_tol) then
               if (check_cap .and. iter >= 1) then
                  if (captors%lands(captor_batch, capture_run, r_curr)) then
                     captured = .true.
                     return
                  end if
               end if
               converged = .true.
               exit
            end if
            ! Flat second order: the third order decides, as in the check
            call lagrangian_third_derivative(lsf, r_curr, lambda_curr, matmul(T_mat, v_min), D3_ddd, &
                                             unreached, eval_error)
            if (allocated(eval_error)) then
               call hand_over(eval_error, error, lsf_error)
               return
            end if
            if (unreached) then
               call fatal_error(error, "Riemannian escape: iterate beyond the level set reach")
               return
            end if
            if (abs(D3_ddd) <= tol_D3) then
               converged = .true.
               exit
            end if
            escaping = .true.
         end if

         if (escaping) then
            ! Unsigned: `grad phi . delta_r` is roundoff at a KKT point
            delta_r = T_mat(:, 1)*v_min(1) + T_mat(:, 2)*v_min(2)
            directional_deriv = 0.0_wp
            alpha = escape_step
         else
            if (mu_min > eig_tol) then
               ! s = -V diag(1/mu) V^T g
               s_tan(1) = -(v_min(1)*g_tan(1) + v_min(2)*g_tan(2))/mu_min*v_min(1) &
                          - (v_max(1)*g_tan(1) + v_max(2)*g_tan(2))/mu_max*v_max(1)
               s_tan(2) = -(v_min(1)*g_tan(1) + v_min(2)*g_tan(2))/mu_min*v_min(2) &
                          - (v_max(1)*g_tan(1) + v_max(2)*g_tan(2))/mu_max*v_max(2)
            else
               s_tan = -g_tan/max(norm_g_tan, 1.0e-10_wp)
            end if
            delta_r = T_mat(:, 1)*s_tan(1) + T_mat(:, 2)*s_tan(2)
            dlen = norm2(delta_r)
            ! Early duplicate capture: an uncapped Newton step on a convex patch
            ! predicts the point it lands on
            if (check_cap .and. iter >= 1 .and. mu_min > eig_tol .and. dlen < max_step) then
               if (captors%lands(captor_batch, capture_run, r_curr + delta_r)) then
                  captured = .true.
                  return
               end if
            end if
            if (dlen > max_step) delta_r = delta_r*(max_step/dlen)
            directional_deriv = dot_product(grad_phi, delta_r)
            alpha = cfg%escape_alpha_init
         end if

         ! An escape step tries both signs and keeps the lower retracted objective;
         ! when `avoid_dir` is given at `iter = 0`, force the opposite 3D sign directly
         nsign = 1
         if (escaping) then
            if (iter == 0 .and. present(avoid_dir)) then
               if (dot_product(delta_r, avoid_dir) > 0.0_wp) delta_r = -delta_r
               nsign = 1
            else
               nsign = 2
            end if
         end if
         ls_success = .false.
         chosen_sign = 1
         have_pos = .false.
         have_neg = .false.
         phi_pos = 0.0_wp
         phi_neg = 0.0_wp
         do ls_iter = 1, escape_max_ls_iter
            have_trial = .false.
            have_pos = .false.
            have_neg = .false.
            do isign = 1, nsign
               r_trial = r_curr + merge(alpha, -alpha, isign == 1)*delta_r
               call retract_to_surface(cfg, phi_alpha, lsf, r_trial, anchor, r_cand, lambda_cand, retract_error, &
                                       lsf_error=eval_error)
               if (allocated(eval_error)) then
                  call hand_over(eval_error, error, lsf_error)
                  return
               end if
               if (allocated(retract_error)) then
                  deallocate (retract_error)
                  cycle
               end if
               phi_cand = 0.5_wp*phi_alpha*sum((r_cand - anchor)**2)
               if (isign == 1) then
                  have_pos = .true.
                  phi_pos = phi_cand
               else
                  have_neg = .true.
                  phi_neg = phi_cand
               end if
               if (have_trial) then
                  if (phi_cand >= phi_trial) cycle
               end if
               r_retracted = r_cand
               lambda_trial = lambda_cand
               phi_trial = phi_cand
               chosen_sign = merge(1, -1, isign == 1)
               have_trial = .true.
            end do
            if (.not. have_trial) then
               alpha = alpha*escape_alpha_reduce
               cycle
            end if

            if (escaping) then
               ! No first-order decrease to demand off a stationary point: the objective must drop
               drop_min = max(16.0_wp*epsilon(phi_curr)*abs(phi_curr), &
                              min(escape_drop_max, 0.25_wp*eig_tol*alpha**2))
               ls_success = phi_trial < phi_curr - drop_min
            else
               ! Sufficient decrease; a bare decrease accepts the overshoot of a long Newton step
               ls_success = phi_trial <= phi_curr + escape_armijo_c*alpha*directional_deriv
            end if
            if (ls_success) exit
            alpha = alpha*escape_alpha_reduce
         end do

         if (iter == 0 .and. escaping .and. ls_success) then
            if (present(escape_dir)) then
               dlen = norm2(delta_r)
               if (dlen > 0.0_wp) escape_dir = (real(chosen_sign, wp)/dlen)*delta_r
            end if
            if (present(is_near_sym)) then
               is_near_sym = have_pos .and. have_neg .and. &
                             abs(phi_pos - phi_neg) <= saddle_sym_rel*max(1.0_wp, phi_init)
            end if
         end if

         if (.not. ls_success) exit

         r_curr = r_retracted
         lambda_curr = lambda_trial
         phi_curr = phi_trial
         r_out = r_curr
         lambda_out = lambda_curr

         if (ctx%debug) then
            call ctx%debug_message("   escape "//to_string(iter)//" "//to_string(phi_curr, "(es12.4)")//" "// &
                                   to_string(norm2(r_curr - anchor), "(f12.6)")//" "// &
                                   to_string(norm_g_tan, "(es12.4)")//" "//to_string(alpha, "(es10.2)")//" "// &
                                   to_string(mu_min, "(es12.4)")//" "//merge("T", "F", escaping))
         end if
      end do

      if (.not. converged .and. iter >= cfg%escape_max_iter) then
         call fatal_error(error, "Riemannian Newton did not converge within max iterations")
      end if
   end subroutine riemannian_newton_escape

   !* ================================================================================= *!
   !*                                 Surface retraction                                *!
   !* ================================================================================= *!

   !> Retract a point onto `S = 0` by Newton steps along the gradient
   !>
   !> - the multiplier is the least-squares stationarity estimate
   !>   `grad phi . grad S / ||grad S||^2` at the retracted point
   !>
   !> @param[in]    cfg        Polish settings
   !> @param[in]    phi_alpha  Anchor objective weight
   !> @param[inout] lsf        Thread-local level set
   !> @param[in]    r_init     Start point
   !> @param[in]    anchor     Anchor
   !> @param[out]   r_out      Point with `|S| < 1e-12`
   !> @param[out]   lambda_out Multiplier at `r_out`
   !> @param[out]   error      Degenerate gradient, unreached trial or iteration budget exhausted
   !> @param[out]   lsf_error  LSF evaluation failure; folded into `error` when absent
   subroutine retract_to_surface(cfg, phi_alpha, lsf, r_init, anchor, r_out, lambda_out, error, lsf_error)
      type(drop_polish_parameters_type), intent(in) :: cfg
      real(wp), intent(in) :: phi_alpha
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: r_init(3)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(out) :: r_out(3)
      real(wp), intent(out) :: lambda_out
      type(error_type), allocatable, intent(out) :: error
      type(error_type), allocatable, intent(out), optional :: lsf_error

      real(wp), parameter :: tol_S = 1.0e-12_wp

      type(error_type), allocatable :: eval_error
      real(wp) :: r_curr(3), S_val, grad_S(3), norm_grad_S2
      integer :: iter
      logical :: unreached

      r_curr = r_init
      r_out = r_init
      lambda_out = 0.0_wp

      do iter = 1, cfg%retract_max_iter
         call point_jet(lsf, r_curr, 1, S_val, unreached, eval_error, grad_S)
         if (allocated(eval_error)) then
            call hand_over(eval_error, error, lsf_error)
            return
         end if
         ! Recoverable: the escape walk backs off its step
         if (unreached) then
            call fatal_error(error, "Retraction failed: trial point beyond the level set reach")
            return
         end if
         norm_grad_S2 = dot_product(grad_S, grad_S)

         if (abs(S_val) < tol_S) then
            if (norm_grad_S2 < 1.0e-14_wp) then
               call fatal_error(error, "Retraction failed: degenerate LSF gradient at converged point")
               return
            end if
            r_out = r_curr
            lambda_out = dot_product(phi_alpha*(r_curr - anchor), grad_S)/norm_grad_S2
            return
         end if

         if (norm_grad_S2 < 1.0e-14_wp) then
            call fatal_error(error, "Retraction failed: degenerate LSF gradient")
            return
         end if
         r_curr = r_curr - (S_val/norm_grad_S2)*grad_S
      end do

      call fatal_error(error, "Retraction to surface did not converge")
   end subroutine retract_to_surface

   !> Third directional derivative of the Lagrangian, `D^3 L[d,d,d] = -lambda S_ijk d_i d_j d_k`
   !>
   !> - `Hess phi` is constant, so only the level set contributes
   !>
   !> @param[inout] lsf       Thread-local level set
   !> @param[in]    r         Surface point
   !> @param[in]    lambda    Multiplier at `r`
   !> @param[in]    d         Tangent direction
   !> @param[out]   D3        Third directional derivative
   !> @param[out]   unreached Point beyond the level set reach; `D3` is then zero and meaningless
   !> @param[out]   error     LSF evaluation failure
   subroutine lagrangian_third_derivative(lsf, r, lambda, d, D3, unreached, error)
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: r(3)
      real(wp), intent(in) :: lambda
      real(wp), intent(in) :: d(3)
      real(wp), intent(out) :: D3
      logical, intent(out) :: unreached
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: S3, grad3(3), hess3(6), third3(10)

      D3 = 0.0_wp
      call point_jet(lsf, r, 3, S3, unreached, error, grad3, hess3, third3)
      if (allocated(error) .or. unreached) return
      D3 = -lambda*packed_full_contract3(third3, d, 3)
   end subroutine lagrangian_third_derivative

   !* ================================================================================= *!
   !*                            Single-point level set jets                            *!
   !* ================================================================================= *!

   !> Level set jet at one point through the batched entry
   !>
   !> - `unreached`: no atom within the screening reach; the outputs are then
   !>   unspecified and must not be read
   !>
   !> @param[inout] lsf       Thread-local level set
   !> @param[in]    point     Evaluation point
   !> @param[in]    order     Highest derivative (0..3)
   !> @param[out]   S         Value
   !> @param[out]   unreached Point beyond the level set reach
   !> @param[out]   error     LSF evaluation failure
   !> @param[out]   grad      Gradient, required for order >= 1
   !> @param[out]   hess      Hessian, packed (6), required for order >= 2
   !> @param[out]   third     Third derivatives, packed (10), required for order 3
   subroutine point_jet(lsf, point, order, S, unreached, error, grad, hess, third)
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: order
      real(wp), intent(out) :: S
      logical, intent(out) :: unreached
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(out), optional :: grad(3)
      real(wp), intent(out), optional :: hess(:)
      real(wp), intent(out), optional :: third(10)

      !> One-column batch; Hessian packed (6, 1), third derivatives packed (10, 1)
      real(wp) :: p(3, 1), v(1), g(3, 1), h(6, 1), t3(10, 1)

      unreached = .false.
      p(:, 1) = point
      select case (order)
      case (0)
         call lsf%evaluate_batch(p, 0, v, error=error)
      case (1)
         call lsf%evaluate_batch(p, 1, v, lsf1_r=g, error=error)
      case (2)
         call lsf%evaluate_batch(p, 2, v, lsf1_r=g, lsf2_rr_packed=h, error=error)
      case default
         call lsf%evaluate_batch(p, 3, v, lsf1_r=g, lsf2_rr_packed=h, lsf3_rrr=t3, error=error)
      end select
      if (allocated(error)) return
      unreached = lsf%unreached(1)
      if (unreached) return
      S = v(1)
      if (order >= 1 .and. present(grad)) grad = g(:, 1)
      if (order >= 2 .and. present(hess)) hess = h(:, 1)
      if (order >= 3 .and. present(third)) third = t3(:, 1)
   end subroutine point_jet

   !> Route an LSF failure to `lsf_error` when the caller asked for it
   !>
   !> @param[inout] eval_error LSF failure, consumed
   !> @param[out]   error      Classification channel
   !> @param[out]   lsf_error  Dedicated LSF channel
   subroutine hand_over(eval_error, error, lsf_error)
      type(error_type), allocatable, intent(inout) :: eval_error
      type(error_type), allocatable, intent(inout) :: error
      type(error_type), allocatable, intent(inout), optional :: lsf_error

      if (present(lsf_error)) then
         call move_alloc(eval_error, lsf_error)
      else
         call move_alloc(eval_error, error)
      end if
   end subroutine hand_over

   !* ================================================================================= *!
   !*                                   Saddle cache                                    *!
   !* ================================================================================= *!

   !> Clear all stored saddles for `1..nkey` while keeping allocated buffers
   !>
   !> @param[inout] self Saddle cache
   !> @param[in]    nkey Number of active anchor keys
   subroutine saddle_cache_reset(self, nkey)
      class(saddle_cache_type), intent(inout) :: self
      integer, intent(in) :: nkey

      integer :: cap

      if (nkey > 0) then
         if (.not. allocated(self%head)) then
            allocate (self%head(max(64, nkey)))
         else if (size(self%head) < nkey) then
            cap = max(nkey, 2*size(self%head))
            deallocate (self%head)
            allocate (self%head(cap))
         end if
         self%head(1:nkey) = 0
      end if
      self%nkey = max(0, nkey)
      self%n = 0
   end subroutine saddle_cache_reset

   !> Find a live saddle of anchor `key` within `sqrt(saddle_dedup_tol2)` of `x`
   !>
   !> @param[in] self Saddle cache
   !> @param[in] key  Anchor key
   !> @param[in] x    Saddle coordinates (3)
   !> @return    idx  Matching saddle index, 0 when none
   pure function saddle_cache_find(self, key, x) result(idx)
      type(saddle_cache_type), intent(in) :: self
      integer, intent(in) :: key
      real(wp), intent(in) :: x(3)
      integer :: idx

      integer :: m
      real(wp) :: d(3)

      idx = 0
      if (key < 1 .or. key > self%nkey) return
      m = self%head(key)
      do while (m > 0)
         d = self%x(:, m) - x
         if (dot_product(d, d) <= saddle_dedup_tol2) then
            idx = m
            return
         end if
         m = self%next(m)
      end do
   end function saddle_cache_find

   !> What a visit of saddle `m` does
   !>
   !> - asymmetric saddle: captured once its first direction is confirmed, else waits
   !> - near-symmetric saddle: walks the opposite sign while no walk took it;
   !>   captured once both directions are confirmed, else waits
   !>
   !> @param[in] self Saddle cache
   !> @param[in] m    Live saddle index
   !> @return    action `visit_capture`, `visit_wait` or `visit_opposite`
   pure function saddle_cache_visit(self, m) result(action)
      type(saddle_cache_type), intent(in) :: self
      integer, intent(in) :: m
      integer :: action

      if (self%near_sym(m) .and. self%lead(2, m) == lead_free) then
         action = visit_opposite
      else if (self%lead(1, m) == lead_confirmed .and. &
               (self%lead(2, m) == lead_confirmed .or. .not. self%near_sym(m))) then
         action = visit_capture
      else
         action = visit_wait
      end if
   end function saddle_cache_visit

   !> Append a newly escaped saddle to anchor `key`'s list, its first direction pending
   !>
   !> @param[inout] self     Saddle cache
   !> @param[in]    key      Anchor key
   !> @param[in]    x        Saddle coordinates (3)
   !> @param[in]    dir      Unit 3D escape direction of the first walk (3)
   !> @param[in]    near_sym Whether the two initial escape signs had near-equal objectives
   !> @param[out]   idx      Index of the new saddle
   subroutine saddle_cache_add(self, key, x, dir, near_sym, idx)
      type(saddle_cache_type), intent(inout) :: self
      integer, intent(in) :: key
      real(wp), intent(in) :: x(3), dir(3)
      logical, intent(in) :: near_sym
      integer, intent(out) :: idx

      integer :: cap

      idx = 0
      if (key < 1) return
      if (key > self%nkey) then
         if (.not. allocated(self%head)) then
            call grow_array(self%head, max(64, key))
         else if (key > size(self%head)) then
            call grow_array(self%head, max(key, 2*size(self%head)))
         end if
         ! Slots past `nkey` may hold lists of an earlier reset
         self%head(self%nkey + 1:key) = 0
         self%nkey = key
      end if

      cap = 0
      if (.not. allocated(self%next)) then
         cap = 64
      else if (self%n + 1 > size(self%next)) then
         cap = 2*size(self%next)
      end if
      if (cap > 0) then
         call grow_array(self%next, cap)
         call grow_array(self%key, cap)
         call grow_array(self%lead, 2, cap)
         call grow_array(self%near_sym, cap)
         call grow_array(self%x, 3, cap)
         call grow_array(self%dir, 3, cap)
      end if

      self%n = self%n + 1
      idx = self%n
      self%next(idx) = self%head(key)
      self%key(idx) = key
      self%x(:, idx) = x
      self%dir(:, idx) = dir
      self%near_sym(idx) = near_sym
      self%lead(:, idx) = [lead_pending, lead_free]
      self%head(key) = idx
   end subroutine saddle_cache_add

   !> Settle the walk leading direction `sign` of saddle `m`
   !>
   !> - confirmed: the direction is `lead_confirmed`
   !> - failed: the saddle leaves its key's list, so a later visit walks as a new
   !>   first visit; a saddle dropped before stays dropped
   !>
   !> @param[inout] self      Saddle cache
   !> @param[in]    m         Saddle index
   !> @param[in]    sign      Direction, 1 first or 2 opposite
   !> @param[in]    confirmed Whether the walk ended accepted or captured by an accepted point
   subroutine saddle_cache_settle(self, m, sign, confirmed)
      type(saddle_cache_type), intent(inout) :: self
      integer, intent(in) :: m
      integer, intent(in) :: sign
      logical, intent(in) :: confirmed

      integer :: p

      if (self%lead(1, m) == lead_dropped) return
      if (confirmed) then
         self%lead(sign, m) = lead_confirmed
         return
      end if
      if (self%head(self%key(m)) == m) then
         self%head(self%key(m)) = self%next(m)
      else
         p = self%head(self%key(m))
         do while (self%next(p) /= m)
            p = self%next(p)
         end do
         self%next(p) = self%next(m)
      end if
      self%lead(:, m) = lead_dropped
   end subroutine saddle_cache_settle

   !> Release saddle cache storage
   !>
   !> @param[inout] self Saddle cache
   subroutine saddle_cache_destroy(self)
      class(saddle_cache_type), intent(inout) :: self

      if (allocated(self%head)) deallocate (self%head)
      if (allocated(self%next)) deallocate (self%next)
      if (allocated(self%key)) deallocate (self%key)
      if (allocated(self%x)) deallocate (self%x)
      if (allocated(self%dir)) deallocate (self%dir)
      if (allocated(self%near_sym)) deallocate (self%near_sym)
      if (allocated(self%lead)) deallocate (self%lead)
      self%nkey = 0
      self%n = 0
   end subroutine saddle_cache_destroy

end module moist_cavity_drop_projection_polish
