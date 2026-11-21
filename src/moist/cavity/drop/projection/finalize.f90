!> Per-anchor finalize of DROP projections: from accepted columns to quadrature rows
!>
!>   1. dedupe within the anchor's run (`branch_weight%sep_cut`)
!>   2. CP Jacobian scaling `cpjac` and gate
!>      `w_f = f_crit * f_foc * f_hor * f_reach` (times the optional
!>      prune switch) per point, from Newton's final jet ([[point_gate]])
!>   3. branches: the unique points in canonical order (coordinates,
!>      `tie_tol = 1e-8`)
!>   4. branch weights over the branches, from their objectives and gates, see
!>      [[moist_cavity_drop_branching]]
!>   5. `wleb = anchor wleb * cpjac`; a branch becomes a row when
!>      `wleb * f_anchor * wbranch > wleb_cut`
!>
!> - the gate and the branch share are one number, `wbranch`, which enters the
!>   row's switching value, `f = f_anchor * wbranch`; `wleb` is the geometric
!>   weight of the patch alone and sets the Gaussian width
!> - the rows of an anchor are the set its weights are formed over: a branch
!>   below the cut is taken out and the weights of the others are formed again,
!>   until every branch left is a row, so the derivatives see every point the
!>   weights depend on and no row has zero area
!> - a point with a closed gate takes no weight and gives none, so leaving it
!>   out is exact; a branch `b` is taken out when
!>   `wleb * cpjac_b * f_anchor * wbranch_b <= wleb_cut`, which can move the
!>   `wbranch` of a sibling `c` by order one but its area by at most
!>   `wleb_cut * cpjac_c / cpjac_b`
!> - a point `branch%dpdist_max` or more behind a fully open one has weight zero
!>   and never is a row, which is what lets the search stop there
!> - `branch` is the position of the row among the rows of its anchor
module moist_cavity_drop_projection_finalize
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_math_linalg, only: setup_tangent_frame
   use moist_math_packed, only: packed_index
   use moist_math_smoothstep, only: moist_math_smoothstep_type
   use moist_cavity_drop_branching, only: branch_weight_type, branch_pdist2_max
   use moist_cavity_drop_projection_batch, only: projection_batch_type, point_accepted
   use moist_cavity_drop_projection_parameters, only: drop_finalize_parameters_type
   implicit none(type, external)
   private

   public :: finalize_anchor_type, projection_rows_type, finalize_scratch_type
   public :: finalize_anchor, dedupe_candidates, canonical_branch_order, cp_jacobian_weight, point_gate

   !> Separation below which two coordinates count as tied in the canonical order
   real(wp), parameter :: tie_tol = 1.0e-8_wp
   !> Determinant below which the tangent KKT block is singular
   real(wp), parameter :: det_B_guard = 1.0e-30_wp
   !> Gradient norm below which a point has no normal
   real(wp), parameter :: g_norm_floor = 1.0e-14_wp


   !> Per-anchor scalars the finalize reads
   type :: finalize_anchor_type
      !> Anchor slot
      integer :: index = 0
      !> Anchor position
      real(wp) :: xyz(3) = 0.0_wp
      !> Owning atom and its center
      integer :: owner = 0
      real(wp) :: owner_center(3) = 0.0_wp
      !> Anchor quadrature weight entering the projection
      real(wp) :: wleb = 0.0_wp
      !> Combined anchor switching value
      real(wp) :: f = 0.0_wp
   end type finalize_anchor_type

   !> Append-only per-thread buffer of quadrature rows, grown geometrically
   type :: projection_rows_type
      !> Rows in use
      integer :: n = 0
      !> Allocated rows
      integer :: capacity = 0
      !> Anchor slot, canonical branch index
      integer, allocatable :: anchor(:), branch(:)
      !> Point and unit normal (3, capacity)
      real(wp), allocatable :: xyz(:, :), normal(:, :)
      !> Multiplier, anchor distance, objective
      real(wp), allocatable :: lambda(:), pdist(:), phi(:)
      !> Geometric quadrature weight; branch weight, the gate included
      real(wp), allocatable :: wleb(:), wbranch(:)
      !> CP Jacobian scaling and gate
      real(wp), allocatable :: cpjac(:), w_f(:)
   contains
      !> Forget all rows, keeping the storage
      procedure :: clear => rows_clear
      !> Grow the storage
      procedure :: reserve => rows_reserve
      !> Release the storage
      procedure :: destroy => rows_destroy
   end type projection_rows_type

   !> Per-thread finalize scratch on the heap; an octree run can hand one anchor
   !> thousands of columns, too many for a worker thread's stack
   type :: finalize_scratch_type
      !> Allocated columns
      integer :: capacity = 0
      !> Accepted, canonical-order and deduplicated column lists
      integer, allocatable :: acc(:), order(:), keep_col(:)
      !> Per unique point: scaling and gate
      real(wp), allocatable :: jac(:), gate(:)
      !> Dedupe keep flags
      logical, allocatable :: keep(:)
      !> Points (3, capacity)
      real(wp), allocatable :: pts(:, :)
      !> Per-branch scaling, gate, objective and weight
      real(wp), allocatable :: cpjac(:), w_f(:), phi(:), wbr(:)
   end type finalize_scratch_type

contains

   !* ================================================================================= *!
   !*                                Per-anchor finalize                                *!
   !* ================================================================================= *!

   !> Turn one anchor's accepted columns into quadrature rows
   !>
   !> - the weights are formed over the rows alone: branches below the cut,
   !>   those with a closed gate first of all, are taken out and the weights
   !>   formed again until none is left below it
   !> - a negative cut keeps every unique point as a row, whatever its weight
   !>
   !> @param[in]    cfg           Finalize settings
   !> @param[in]    phi_alpha     Anchor objective weight
   !> @param[in]    batch         Batch holding the anchor's columns
   !> @param[in]    first         First column of the anchor's run
   !> @param[in]    last          Last column of the anchor's run; non-accepted columns are skipped
   !> @param[in]    anchor        Per-anchor scalars
   !> @param[in]    f_crit        Gradient-norm switch
   !> @param[in]    f_foc         Focal switch on the focal margin `beta_min(B)/phi_alpha`
   !> @param[in]    f_hor         Horizon switch on the signed `n_surf . n_sph`
   !> @param[in]    f_wleb        Prune switch, read when `cfg%wleb_prune`
   !> @param[in]    f_reach       Reach gate on the projection distance
   !> @param[in]    branch_weight Branch weight model and branch separation
   !> @param[inout] rows          Per-thread row buffer
   !> @param[out]   npoint        Unique accepted points (0: the anchor is lost)
   !> @param[out]   nrows         Rows appended
   !> @param[inout] sc            Per-thread scratch
   !> @param[out]   error         Singular tangent Jacobian (build-fatal)
   subroutine finalize_anchor(cfg, phi_alpha, batch, first, last, anchor, f_crit, f_foc, f_hor, f_wleb, &
                              f_reach, branch_weight, rows, npoint, nrows, sc, error)
      type(drop_finalize_parameters_type), intent(in) :: cfg
      real(wp), intent(in) :: phi_alpha
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: first, last
      type(finalize_anchor_type), intent(in) :: anchor
      class(moist_math_smoothstep_type), intent(in) :: f_crit
      class(moist_math_smoothstep_type), intent(in) :: f_foc
      class(moist_math_smoothstep_type), intent(in) :: f_hor
      class(moist_math_smoothstep_type), intent(in) :: f_wleb
      class(moist_math_smoothstep_type), intent(in) :: f_reach
      type(branch_weight_type), intent(in) :: branch_weight
      type(projection_rows_type), intent(inout) :: rows
      integer, intent(out) :: npoint
      integer, intent(out) :: nrows
      type(finalize_scratch_type), intent(inout) :: sc
      type(error_type), allocatable, intent(out) :: error

      integer :: nbranch, nkeep, u, b, j, row

      nrows = 0
      call gather_points(cfg, phi_alpha, batch, first, last, anchor, f_crit, f_foc, f_hor, f_wleb, f_reach, &
                         branch_weight%sep_cut, npoint, sc, error)
      if (allocated(error)) return
      if (npoint == 0) return

      nbranch = npoint
      do u = 1, npoint
         sc%acc(u) = u
         sc%pts(:, u) = batch%x(:, sc%keep_col(u))
      end do

      !> Canonical order; `acc` then holds the batch columns
      call canonical_branch_order(sc%pts(:, 1:nbranch), sc%order(1:nbranch))
      do b = 1, nbranch
         u = sc%acc(sc%order(b))
         sc%phi(b) = batch%phi(sc%keep_col(u))
         sc%cpjac(b) = sc%jac(u)
         sc%w_f(b) = sc%gate(u)
         sc%order(b) = sc%keep_col(u)
      end do

      !> Weights over the branches above the cut: take the ones below it out
      !> until the set is stable; removal keeps the canonical order
      do
         call branch_weight%weights(sc%phi(1:nbranch), sc%w_f(1:nbranch), sc%wbr(1:nbranch))
         nkeep = 0
         do b = 1, nbranch
            if (.not. anchor%wleb*sc%cpjac(b)*anchor%f*sc%wbr(b) > cfg%wleb_cut) cycle
            nkeep = nkeep + 1
            if (nkeep == b) cycle
            sc%phi(nkeep) = sc%phi(b)
            sc%cpjac(nkeep) = sc%cpjac(b)
            sc%w_f(nkeep) = sc%w_f(b)
            sc%wbr(nkeep) = sc%wbr(b)
            sc%order(nkeep) = sc%order(b)
         end do
         if (nkeep == nbranch) exit
         nbranch = nkeep
         if (nbranch == 0) return
      end do

      do b = 1, nbranch
         j = sc%order(b)
         if (rows%n + 1 > rows%capacity) call rows%reserve(rows%n + 1)
         rows%n = rows%n + 1
         row = rows%n
         rows%anchor(row) = anchor%index
         rows%branch(row) = b
         rows%xyz(:, row) = batch%x(:, j)
         rows%normal(:, row) = batch%normal(:, j)
         rows%lambda(row) = batch%lambda(j)
         rows%pdist(row) = batch%pdist(j)
         rows%phi(row) = sc%phi(b)
         rows%wleb(row) = anchor%wleb*sc%cpjac(b)
         rows%wbranch(row) = sc%wbr(b)
         rows%cpjac(row) = sc%cpjac(b)
         rows%w_f(row) = sc%w_f(b)
         nrows = nrows + 1
      end do
   end subroutine finalize_anchor

   !> Unique accepted points of one anchor with their scalings and gates
   !>
   !> - on return `sc%keep_col(1:npoint)` holds the columns in run order and
   !>   `sc%jac`, `sc%gate` their values
   !>
   !> @param[in]    cfg       Finalize settings
   !> @param[in]    phi_alpha Anchor objective weight
   !> @param[in]    batch     Batch holding the anchor's columns
   !> @param[in]    first     First column of the anchor's run
   !> @param[in]    last      Last column of the anchor's run
   !> @param[in]    anchor    Per-anchor scalars
   !> @param[in]    f_crit    Gradient-norm switch
   !> @param[in]    f_foc     Focal switch
   !> @param[in]    f_hor     Horizon switch
   !> @param[in]    f_wleb    Prune switch
   !> @param[in]    f_reach   Reach gate
   !> @param[in]    sep_cut   Separation below which two points are one branch (Bohr)
   !> @param[out]   npoint    Unique accepted points
   !> @param[inout] sc        Per-thread scratch
   !> @param[out]   error     Singular tangent Jacobian (build-fatal)
   subroutine gather_points(cfg, phi_alpha, batch, first, last, anchor, f_crit, f_foc, f_hor, f_wleb, f_reach, &
                            sep_cut, npoint, sc, error)
      type(drop_finalize_parameters_type), intent(in) :: cfg
      real(wp), intent(in) :: phi_alpha
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: first, last
      type(finalize_anchor_type), intent(in) :: anchor
      class(moist_math_smoothstep_type), intent(in) :: f_crit
      class(moist_math_smoothstep_type), intent(in) :: f_foc
      class(moist_math_smoothstep_type), intent(in) :: f_hor
      class(moist_math_smoothstep_type), intent(in) :: f_wleb
      class(moist_math_smoothstep_type), intent(in) :: f_reach
      real(wp), intent(in) :: sep_cut
      integer, intent(out) :: npoint
      type(finalize_scratch_type), intent(inout) :: sc
      type(error_type), allocatable, intent(out) :: error

      integer :: nacc, k, u

      npoint = 0
      call reserve_finalize_scratch(sc, last - first + 1)

      !> Accepted columns, then dedupe in run order
      nacc = 0
      do k = first, last
         if (batch%status(k) /= point_accepted) cycle
         nacc = nacc + 1
         sc%acc(nacc) = k
         sc%pts(:, nacc) = batch%x(:, k)
      end do
      if (nacc == 0) return
      call unique_points(sc%pts(:, 1:nacc), sep_cut, sc%keep(1:nacc))
      do k = 1, nacc
         if (.not. sc%keep(k)) cycle
         npoint = npoint + 1
         sc%keep_col(npoint) = sc%acc(k)
      end do

      do u = 1, npoint
         call point_gate(cfg, phi_alpha, batch, sc%keep_col(u), anchor, f_crit, f_foc, f_hor, f_wleb, f_reach, &
                         sc%jac(u), sc%gate(u), error)
         if (allocated(error)) return
      end do
   end subroutine gather_points

   !> CP Jacobian scaling and full gate of one solved column
   !>
   !> The one evaluation of the gate of a projected point: the search reads it
   !> for a just-accepted point to decide how far it must look, the finalize
   !> for the weights
   !>
   !> - `gate = f_crit * f_foc * f_hor * f_reach`, times the prune
   !>   switch of the geometric weight when `cfg%wleb_prune`
   !> - an anchor below the switching cut has every gate closed
   !>
   !> @param[in]  cfg       Finalize settings
   !> @param[in]  phi_alpha Anchor objective weight
   !> @param[in]  batch     Batch holding the column
   !> @param[in]  j         Solved column
   !> @param[in]  anchor    Per-anchor scalars
   !> @param[in]  f_crit    Gradient-norm switch
   !> @param[in]  f_foc     Focal switch
   !> @param[in]  f_hor     Horizon switch
   !> @param[in]  f_wleb    Prune switch
   !> @param[in]  f_reach   Reach gate on the projection distance
   !> @param[out] jac       Area scaling, 1 under a closed gate
   !> @param[out] gate      Gate
   !> @param[out] error     Singular tangent Jacobian (build-fatal)
   subroutine point_gate(cfg, phi_alpha, batch, j, anchor, f_crit, f_foc, f_hor, f_wleb, f_reach, jac, gate, error)
      type(drop_finalize_parameters_type), intent(in) :: cfg
      real(wp), intent(in) :: phi_alpha
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: j
      type(finalize_anchor_type), intent(in) :: anchor
      class(moist_math_smoothstep_type), intent(in) :: f_crit
      class(moist_math_smoothstep_type), intent(in) :: f_foc
      class(moist_math_smoothstep_type), intent(in) :: f_hor
      class(moist_math_smoothstep_type), intent(in) :: f_wleb
      class(moist_math_smoothstep_type), intent(in) :: f_reach
      real(wp), intent(out) :: jac
      real(wp), intent(out) :: gate
      type(error_type), allocatable, intent(out) :: error

      !> Path from the anchor to the point
      real(wp) :: v(3)

      jac = 1.0_wp
      gate = 0.0_wp
      if (anchor%f < cfg%wleb_cut) return
      call cp_jacobian_weight(phi_alpha, batch%gradS(:, j), batch%hessS(:, j), &
                              batch%lambda(j), anchor%xyz, anchor%owner_center, f_crit, f_foc, &
                              f_hor, jac, gate, error)
      if (allocated(error)) return
      if (gate > 0.0_wp) then
         v = batch%x(:, j) - anchor%xyz
         gate = gate*f_reach%f0(sqrt(dot_product(v, v)))
      end if
      ! The prune switch of the geometric weight is part of the gate, so a
      ! pruned branch gives its share back before it is cut
      if (cfg%wleb_prune) gate = gate*f_wleb%f0(abs(anchor%wleb*jac))
   end subroutine point_gate

   !> Grow the finalize scratch to `n` columns
   !>
   !> @param[inout] sc Scratch
   !> @param[in]    n  Columns to hold
   subroutine reserve_finalize_scratch(sc, n)
      type(finalize_scratch_type), intent(inout) :: sc
      integer, intent(in) :: n

      integer :: cap

      if (n <= sc%capacity) return
      cap = max(16, 2*sc%capacity, n)
      if (allocated(sc%acc)) then
         deallocate (sc%acc, sc%order, sc%keep_col, sc%keep, sc%pts, sc%cpjac, &
                     sc%w_f, sc%phi, sc%wbr, sc%jac, sc%gate)
      end if
      allocate (sc%acc(cap), sc%order(cap), sc%keep_col(cap), sc%keep(cap))
      allocate (sc%pts(3, cap), sc%cpjac(cap), sc%w_f(cap), sc%phi(cap), sc%wbr(cap))
      allocate (sc%jac(cap), sc%gate(cap))
      sc%capacity = cap
   end subroutine reserve_finalize_scratch

   !* ================================================================================= *!
   !*                         Deduplication and canonical order                         *!
   !* ================================================================================= *!

   !> Keep flags for unique points inside the squared-distance window
   !>
   !> - a point within `sep_cut` of an earlier kept point is a duplicate
   !> - a point beyond the window bound over the closest kept one is dropped,
   !>   see [[branch_pdist2_max]]
   !>
   !> @param[in]  pts          Points (3, n)
   !> @param[in]  anchor       Anchor
   !> @param[in]  sep_cut      Duplicate separation (Bohr)
   !> @param[in]  pdist2_slack Constant squared-distance part of the window (Bohr^2)
   !> @param[out] keep         Keep flags (n)
   !> @param[in]  pdist_slack  Optional distance part of the window (Bohr)
   !> @param[in]  pdist_soft   Optional softening length of the distance part (Bohr)
   pure subroutine dedupe_candidates(pts, anchor, sep_cut, pdist2_slack, keep, pdist_slack, pdist_soft)
      real(wp), intent(in) :: pts(:, :)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: sep_cut
      real(wp), intent(in) :: pdist2_slack
      logical, intent(out) :: keep(:)
      real(wp), intent(in), optional :: pdist_slack
      real(wp), intent(in), optional :: pdist_soft

      integer :: i
      real(wp) :: pdist2_min, pdist2_max, d2, dpdist, soft

      call unique_points(pts, sep_cut, keep)

      pdist2_min = huge(1.0_wp)
      do i = 1, size(pts, 2)
         if (keep(i)) pdist2_min = min(pdist2_min, sum((pts(:, i) - anchor)**2))
      end do
      if (pdist2_min == huge(1.0_wp)) return
      dpdist = 0.0_wp
      soft = 0.0_wp
      if (present(pdist_slack)) dpdist = pdist_slack
      if (present(pdist_soft)) soft = pdist_soft
      pdist2_max = branch_pdist2_max(pdist2_min, pdist2_slack, dpdist, soft)
      do i = 1, size(pts, 2)
         if (.not. keep(i)) cycle
         d2 = sum((pts(:, i) - anchor)**2)
         if (d2 > pdist2_max) keep(i) = .false.
      end do
   end subroutine dedupe_candidates

   !> Keep flags for unique points
   !>
   !> - a point within `sep_cut` of an earlier kept point is a duplicate
   !>
   !> @param[in]  pts     Points (3, n)
   !> @param[in]  sep_cut Duplicate separation (Bohr)
   !> @param[out] keep    Keep flags (n)
   pure subroutine unique_points(pts, sep_cut, keep)
      real(wp), intent(in) :: pts(:, :)
      real(wp), intent(in) :: sep_cut
      logical, intent(out) :: keep(:)

      integer :: i, j, m, n_kept
      ! Kept points so far, in run order
      integer, allocatable :: kept(:)

      allocate (kept(size(pts, 2)))
      n_kept = 0
      keep = .false.
      do i = 1, size(pts, 2)
         keep(i) = .true.
         do m = 1, n_kept
            j = kept(m)
            if (sum((pts(:, i) - pts(:, j))**2) < sep_cut**2) then
               keep(i) = .false.
               exit
            end if
         end do
         if (.not. keep(i)) cycle
         n_kept = n_kept + 1
         kept(n_kept) = i
      end do
   end subroutine unique_points

   !> Geometry-derived branch order, so the branch index names the same surface
   !> point on every build, whatever the input order
   !>
   !> - lexicographic on x, y, z, a coordinate compared through its tie cluster:
   !>   values chained by gaps of at most `tie_tol` share one, and the next
   !>   coordinate decides inside it instead of rounding noise
   !> - y clusters are formed inside an x cluster; exact z, x, y break the rest
   !> - clusters, the transitive closure of the pairwise tie, keep the order a
   !>   strict weak one; the pairwise tie alone is not transitive along a chain
   !>   of near-ties and made the order depend on the input order
   !>
   !> @param[in]  pts   Points (3, n)
   !> @param[out] order Permutation into canonical order
   pure subroutine canonical_branch_order(pts, order)
      real(wp), intent(in) :: pts(:, :)
      integer, intent(out) :: order(:)

      !> Tie cluster of x; of y inside the x cluster, numbered across clusters
      integer, allocatable :: cx(:), cy(:)
      integer :: i, n

      n = size(pts, 2)
      do i = 1, n
         order(i) = i
      end do
      if (n < 2) return
      allocate (cx(n), cy(n))

      call sort_by(order, 1)
      cx(order(1)) = 1
      do i = 2, n
         cx(order(i)) = cx(order(i - 1))
         if (pts(1, order(i)) - pts(1, order(i - 1)) > tie_tol) cx(order(i)) = cx(order(i)) + 1
      end do

      call sort_by(order, 2)
      cy(order(1)) = 1
      do i = 2, n
         cy(order(i)) = cy(order(i - 1))
         if (cx(order(i)) /= cx(order(i - 1)) .or. pts(2, order(i)) - pts(2, order(i - 1)) > tie_tol) then
            cy(order(i)) = cy(order(i)) + 1
         end if
      end do

      call sort_by(order, 3)

   contains

      !> Stable insertion sort of `order` by the key of a stage
      pure subroutine sort_by(order, stage)
         integer, intent(inout) :: order(:)
         integer, intent(in) :: stage
         integer :: i, j, moved

         do i = 2, size(order)
            moved = order(i)
            j = i - 1
            do while (j >= 1)
               if (.not. precedes(moved, order(j), stage)) exit
               order(j + 1) = order(j)
               j = j - 1
            end do
            order(j + 1) = moved
         end do
      end subroutine sort_by

      !> Whether point `a` sorts before point `b`: stage 1 on x, stage 2 on
      !> (x cluster, y), stage 3 on (y cluster, z, x, y)
      pure logical function precedes(a, b, stage) result(before)
         integer, intent(in) :: a, b, stage
         !> Exact tie-break order of stage 3: z, x, y
         integer, parameter :: tail(3) = [3, 1, 2]
         integer :: k

         select case (stage)
         case (1)
            before = pts(1, a) < pts(1, b)
         case (2)
            before = cx(a) < cx(b) .or. (cx(a) == cx(b) .and. pts(2, a) < pts(2, b))
         case default
            before = cy(a) < cy(b)
            if (cy(a) /= cy(b)) return
            do k = 1, 3
               before = pts(tail(k), a) < pts(tail(k), b)
               if (pts(tail(k), a) /= pts(tail(k), b)) return
            end do
         end select
      end function precedes

   end subroutine canonical_branch_order

   !* ================================================================================= *!
   !*                           Closest-point Jacobian scaling                          *!
   !* ================================================================================= *!

   !> Area scaling of the anchor-to-surface map and the local gate of the point
   !>
   !>     B = Q^T (w I - lambda Hess S) Q,  tau_k = Q^T t_k,  y_k = w Q B^-1 tau_k,  J = |y1 x y2|
   !>
   !> - `w_f = f_crit(||grad S||) * f_foc(beta_min(B)/w) * f_hor(n_surf . n_sph)`;
   !>   a closed gate leaves `J = 1`
   !> - `beta_min(B)/w` is the focal margin, one on a flat surface and zero at a
   !>   focal point, whatever the objective weight
   !> - `f_hor` closes the gate of a point whose surface normal turns tangent to
   !>   the anchor sphere, before `J` reaches zero, or points into it
   !>
   !> @param[in]  w            Objective weight
   !> @param[in]  g            Level set gradient at the point
   !> @param[in]  H            Level set Hessian at the point, packed (6)
   !> @param[in]  lambda       Multiplier
   !> @param[in]  anchor       Anchor
   !> @param[in]  owner_center Owning atom center
   !> @param[in]  f_crit       Gradient-norm switch
   !> @param[in]  f_foc        Focal switch
   !> @param[in]  f_hor        Horizon switch
   !> @param[out] J            Area scaling
   !> @param[out] w_f          Gate without the path gates and the prune switch
   !> @param[out] error        Singular tangent block after switching
   subroutine cp_jacobian_weight(w, g, H, lambda, anchor, owner_center, f_crit, f_foc, f_hor, J, w_f, error)
      real(wp), intent(in) :: w, g(3), H(:), lambda, anchor(3), owner_center(3)
      class(moist_math_smoothstep_type), intent(in) :: f_crit
      class(moist_math_smoothstep_type), intent(in) :: f_foc
      class(moist_math_smoothstep_type), intent(in) :: f_hor
      real(wp), intent(out) :: J
      real(wp), intent(out) :: w_f
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: A(3, 3), g_norm, n_surf(3), q1(3), q2(3), B11, B12, B22
      real(wp) :: tr_B, det_B, sqrt_disc, beta2, Binv11, Binv12, Binv22
      real(wp) :: n_sph(3), t1(3), t2(3), tau1(2), tau2(2), w1(2), w2(2), y1(3), y2(3), c(3)
      integer :: i, k

      J = 1.0_wp
      g_norm = norm2(g)
      w_f = f_crit%f0(g_norm)
      if (g_norm <= g_norm_floor) return

      do k = 1, 3
         do i = 1, 3
            A(i, k) = -lambda*H(packed_index(i, k))
         end do
      end do
      do i = 1, 3
         A(i, i) = A(i, i) + w
      end do

      n_surf = g/g_norm
      call setup_tangent_frame(n_surf, q1, q2)
      B11 = dot_product(q1, matmul(A, q1))
      B12 = dot_product(q1, matmul(A, q2))
      B22 = dot_product(q2, matmul(A, q2))

      ! Sum of squares rather than tr^2/4 - det, which cancels at umbilics
      tr_B = B11 + B22
      det_B = B11*B22 - B12*B12
      sqrt_disc = hypot(0.5_wp*(B11 - B22), B12)
      beta2 = 0.5_wp*tr_B - sqrt_disc

      ! The switches reach exactly zero; an open gate implies a regular tangent block
      w_f = w_f*f_foc%f0(beta2/w)
      if (w_f == 0.0_wp) return

      n_sph = anchor - owner_center
      w_f = w_f*f_hor%f0(dot_product(n_surf, n_sph)/norm2(n_sph))
      if (w_f == 0.0_wp) return

      if (abs(det_B) <= det_B_guard) then
         call fatal_error(error, "Tangent Jacobian matrix B is singular after switching")
         return
      end if
      Binv11 = B22/det_B
      Binv12 = -B12/det_B
      Binv22 = B11/det_B

      call setup_tangent_frame(n_sph, t1, t2)
      tau1 = [dot_product(q1, t1), dot_product(q2, t1)]
      tau2 = [dot_product(q1, t2), dot_product(q2, t2)]
      w1 = [Binv11*tau1(1) + Binv12*tau1(2), Binv12*tau1(1) + Binv22*tau1(2)]
      w2 = [Binv11*tau2(1) + Binv12*tau2(2), Binv12*tau2(1) + Binv22*tau2(2)]
      y1 = w*(w1(1)*q1 + w1(2)*q2)
      y2 = w*(w2(1)*q1 + w2(2)*q2)
      c(1) = y1(2)*y2(3) - y1(3)*y2(2)
      c(2) = y1(3)*y2(1) - y1(1)*y2(3)
      c(3) = y1(1)*y2(2) - y1(2)*y2(1)
      J = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
   end subroutine cp_jacobian_weight

   !* ================================================================================= *!
   !*                                     Row buffer                                    *!
   !* ================================================================================= *!

   !> Forget all rows, keeping the storage
   !>
   !> @param[inout] self Row buffer
   subroutine rows_clear(self)
      class(projection_rows_type), intent(inout) :: self

      self%n = 0
   end subroutine rows_clear

   !> Grow the storage to at least `required` rows, keeping `1:n`
   !>
   !> @param[inout] self     Row buffer
   !> @param[in]    required Minimum capacity
   subroutine rows_reserve(self, required)
      class(projection_rows_type), intent(inout) :: self
      integer, intent(in) :: required

      integer :: cap

      if (required <= self%capacity) return
      cap = max(16, 2*self%capacity, required)
      call grow_i(self%anchor)
      call grow_i(self%branch)
      call grow_r2(self%xyz)
      call grow_r2(self%normal)
      call grow_r(self%lambda)
      call grow_r(self%pdist)
      call grow_r(self%phi)
      call grow_r(self%wleb)
      call grow_r(self%wbranch)
      call grow_r(self%cpjac)
      call grow_r(self%w_f)
      self%capacity = cap

   contains

      !> Grow one real vector
      subroutine grow_r(a)
         real(wp), allocatable, intent(inout) :: a(:)
         real(wp), allocatable :: tmp(:)

         allocate (tmp(cap))
         if (allocated(a) .and. self%n > 0) tmp(1:self%n) = a(1:self%n)
         call move_alloc(tmp, a)
      end subroutine grow_r

      !> Grow one (3, n) real array
      subroutine grow_r2(a)
         real(wp), allocatable, intent(inout) :: a(:, :)
         real(wp), allocatable :: tmp(:, :)

         allocate (tmp(3, cap))
         if (allocated(a) .and. self%n > 0) tmp(:, 1:self%n) = a(:, 1:self%n)
         call move_alloc(tmp, a)
      end subroutine grow_r2

      !> Grow one integer vector
      subroutine grow_i(a)
         integer, allocatable, intent(inout) :: a(:)
         integer, allocatable :: tmp(:)

         allocate (tmp(cap))
         if (allocated(a) .and. self%n > 0) tmp(1:self%n) = a(1:self%n)
         call move_alloc(tmp, a)
      end subroutine grow_i

   end subroutine rows_reserve

   !> Release the storage
   !>
   !> @param[inout] self Row buffer
   subroutine rows_destroy(self)
      class(projection_rows_type), intent(inout) :: self

      if (allocated(self%anchor)) deallocate (self%anchor)
      if (allocated(self%branch)) deallocate (self%branch)
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%normal)) deallocate (self%normal)
      if (allocated(self%lambda)) deallocate (self%lambda)
      if (allocated(self%pdist)) deallocate (self%pdist)
      if (allocated(self%phi)) deallocate (self%phi)
      if (allocated(self%wleb)) deallocate (self%wleb)
      if (allocated(self%wbranch)) deallocate (self%wbranch)
      if (allocated(self%cpjac)) deallocate (self%cpjac)
      if (allocated(self%w_f)) deallocate (self%w_f)
      self%n = 0
      self%capacity = 0
   end subroutine rows_destroy

end module moist_cavity_drop_projection_finalize
