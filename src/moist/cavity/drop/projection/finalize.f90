!> Per-anchor finalize of DROP projections: from accepted columns to quadrature rows
!>
!> Same order as the old post-processing chain, so areas carry over:
!>   1. dedupe within the anchor's run (`branch_sep_cut`, then the rho^2 window)
!>   2. canonical branch order (coordinates, `tie_tol = 1e-8`)
!>   3. CP Jacobian scaling and `w_f0` per branch, from Newton's final jet
!>   4. `wleb = anchor wleb * cpjac * w_f0`, then the optional prune switch
!>   5. softmax branch weights over *all* deduplicated branches, sibling cut at
!>      `wleb_cut`, re-softmax over the kept siblings, `wleb *= wbranch`
!>   6. rows with `wleb * f <= wleb_cut` are never appended
!>
!> - `branch` keeps the canonical position from step 2, so a row's branch
!>   index can skip the siblings the final cut removed
module moist_cavity_drop_projection_finalize
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_math_linalg, only: setup_tangent_frame
   use moist_math_packed, only: packed_index
   use moist_math_smoothstep, only: moist_math_smoothstep_type
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_projection_batch, only: projection_batch_type, point_accepted
   implicit none(type, external)
   private

   public :: finalize_config_type, finalize_anchor_type, projection_rows_type, finalize_scratch_type
   public :: finalize_anchor, dedupe_candidates, canonical_branch_order, cp_jacobian_weight
   public :: branch_weights_with_cut

   !> Separation below which two coordinates count as tied in the canonical order
   real(wp), parameter :: tie_tol = 1.0e-8_wp
   !> Determinant below which the tangent KKT block is singular
   real(wp), parameter :: det_B_guard = 1.0e-30_wp

   !> Finalize settings
   ! TODO: good idea to introduce a new type here with drop parameters again?
   type :: finalize_config_type
      !> Anchor objective weight
      real(wp) :: phi_alpha = 0.5_wp
      !> Separation below which two feet are one branch (Bohr)
      real(wp) :: branch_sep_cut = 1.0e-8_wp
      !> Squared-distance window over the closest branch (Bohr^2)
      real(wp) :: rho2_slack = 0.0_wp
      !> Quadrature weight cutoff
      real(wp) :: wleb_cut = 5.0e-12_wp
      !> Prune-switch level; 0 disables it
      integer :: wleb_prune_level = 0
   end type finalize_config_type

   !> Per-anchor scalars the finalize reads
   type :: finalize_anchor_type
      !> Anchor slot
      integer :: index = 0
      !> Anchor position
      real(wp) :: xyz(3) = 0.0_wp
      !> Owning atom center
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
      !> Foot and unit normal (3, capacity)
      real(wp), allocatable :: xyz(:, :), normal(:, :)
      !> Multiplier, anchor distance, objective
      real(wp), allocatable :: lambda(:), rho(:), phi(:)
      !> Final quadrature weight, branch weight
      real(wp), allocatable :: wleb(:), wbranch(:)
      !> CP Jacobian scaling and switching weight
      real(wp), allocatable :: cpjac(:), w_f(:)
      !> Distance to the owning atom
      real(wp), allocatable :: r_iI(:)
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
      !> Dedupe and sibling-cut keep flags
      logical, allocatable :: keep(:), keep_w(:)
      !> Feet (3, capacity)
      real(wp), allocatable :: pts(:, :)
      !> Per-branch weight, scaling, switching, objective, branch weight
      real(wp), allocatable :: wleb(:), cpjac(:), w_f(:), phi(:), wbr(:)
      !> Sibling-cut workspace (2 * capacity)
      real(wp), allocatable :: wkept(:)
   end type finalize_scratch_type

contains

   !* ================================================================================= *!
   !*                                Per-anchor finalize                                *!
   !* ================================================================================= *!

   !> Turn one anchor's accepted columns into quadrature rows
   !>
   !> @param[in]    cfg           Finalize settings
   !> @param[in]    batch         Batch holding the anchor's columns
   !> @param[in]    first         First column of the anchor's run
   !> @param[in]    last          Last column of the anchor's run; non-accepted columns are skipped
   !> @param[in]    anchor        Per-anchor scalars
   !> @param[in]    f_crit        Gradient-norm switch
   !> @param[in]    f_foc         Focal switch on the smallest tangent eigenvalue
   !> @param[in]    f_wleb        Prune switch, read when `wleb_prune_level > 0`
   !> @param[in]    branch_weight Softmax branch weight model
   !> @param[inout] rows          Per-thread row buffer
   !> @param[out]   nbranch       Branches after dedupe (0: the anchor is lost)
   !> @param[out]   nrows         Rows appended
   !> @param[inout] sc            Per-thread scratch
   !> @param[out]   error         Singular tangent Jacobian (build-fatal)
   subroutine finalize_anchor(cfg, batch, first, last, anchor, f_crit, f_foc, f_wleb, branch_weight, &
                              rows, nbranch, nrows, sc, error)
      type(finalize_config_type), intent(in) :: cfg
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: first, last
      type(finalize_anchor_type), intent(in) :: anchor
      class(moist_math_smoothstep_type), intent(in) :: f_crit
      class(moist_math_smoothstep_type), intent(in) :: f_foc
      class(moist_math_smoothstep_type), intent(in) :: f_wleb
      type(branch_weight_type), intent(in) :: branch_weight
      type(projection_rows_type), intent(inout) :: rows
      integer, intent(out) :: nbranch
      integer, intent(out) :: nrows
      type(finalize_scratch_type), intent(inout) :: sc
      type(error_type), allocatable, intent(out) :: error

      integer :: nacc, k, b, j, n_kept, row

      nbranch = 0
      nrows = 0
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
      call dedupe_candidates(sc%pts(:, 1:nacc), anchor%xyz, cfg%branch_sep_cut, cfg%rho2_slack, sc%keep(1:nacc))
      do k = 1, nacc
         if (.not. sc%keep(k)) cycle
         nbranch = nbranch + 1
         sc%keep_col(nbranch) = sc%acc(k)
         sc%pts(:, nbranch) = batch%x(:, sc%acc(k))
      end do
      if (nbranch == 0) return

      !> Canonical order
      call canonical_branch_order(sc%pts(:, 1:nbranch), sc%order(1:nbranch))

      !> CP Jacobian, switching and prune per branch
      do b = 1, nbranch
         j = sc%keep_col(sc%order(b))
         sc%phi(b) = batch%phi(j)
         if (anchor%f < cfg%wleb_cut) then
            sc%cpjac(b) = 1.0_wp
            sc%w_f(b) = 0.0_wp
         else
            call cp_jacobian_weight(cfg%phi_alpha, batch%gradS(:, j), batch%hessS(:, j), &
                                    batch%lambda(j), anchor%xyz, anchor%owner_center, f_crit, f_foc, &
                                    cfg%wleb_cut, sc%cpjac(b), sc%w_f(b), error)
            if (allocated(error)) return
         end if
         sc%wleb(b) = anchor%wleb*sc%cpjac(b)*sc%w_f(b)
         if (cfg%wleb_prune_level > 0) sc%wleb(b) = sc%wleb(b)*f_wleb%f0(abs(sc%wleb(b)))
      end do

      !> Softmax branch weights with the sibling cut
      call branch_weights_with_cut(branch_weight, sc%phi(1:nbranch), cfg%wleb_cut, sc%wbr(1:nbranch), &
                                   n_kept, sc%wkept(1:2*nbranch), sc%keep_w(1:nbranch))
      if (nbranch > 1) sc%wleb(1:nbranch) = sc%wleb(1:nbranch)*sc%wbr(1:nbranch)

      !> Cut and append
      do b = 1, nbranch
         if (.not. sc%wleb(b)*anchor%f > cfg%wleb_cut) cycle
         j = sc%keep_col(sc%order(b))
         if (rows%n + 1 > rows%capacity) call rows%reserve(rows%n + 1)
         rows%n = rows%n + 1
         row = rows%n
         rows%anchor(row) = anchor%index
         rows%branch(row) = b
         rows%xyz(:, row) = batch%x(:, j)
         rows%normal(:, row) = batch%normal(:, j)
         rows%lambda(row) = batch%lambda(j)
         rows%rho(row) = batch%rho(j)
         rows%phi(row) = sc%phi(b)
         rows%wleb(row) = sc%wleb(b)
         rows%wbranch(row) = sc%wbr(b)
         rows%cpjac(row) = sc%cpjac(b)
         rows%w_f(row) = sc%w_f(b)
         rows%r_iI(row) = norm2(batch%x(:, j) - anchor%owner_center)
         nrows = nrows + 1
      end do
   end subroutine finalize_anchor

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
         deallocate (sc%acc, sc%order, sc%keep_col, sc%keep, sc%keep_w, sc%pts, sc%wleb, sc%cpjac, &
                     sc%w_f, sc%phi, sc%wbr, sc%wkept)
      end if
      allocate (sc%acc(cap), sc%order(cap), sc%keep_col(cap), sc%keep(cap), sc%keep_w(cap))
      allocate (sc%pts(3, cap), sc%wleb(cap), sc%cpjac(cap), sc%w_f(cap), sc%phi(cap), sc%wbr(cap))
      allocate (sc%wkept(2*cap))
      sc%capacity = cap
   end subroutine reserve_finalize_scratch

   !> Softmax branch weights over all siblings with the sibling cut
   !>
   !> - one sibling: weight 1
   !> - siblings whose softmax weight is not above `wleb_cut` get weight 0;
   !>   when that would drop all of them the strongest one is kept
   !> - the kept siblings are re-softmaxed so their weights sum to one
   !>
   !> @param[in]  branch_weight Softmax model
   !> @param[in]  phi           Objective per sibling
   !> @param[in]  wleb_cut      Weight cutoff
   !> @param[out] wbr           Branch weight per sibling, 0 when cut
   !> @param[out] n_kept        Siblings kept
   !> @param[out] work          Workspace, at least `2 size(phi)`
   !> @param[out] keep          Workspace, at least `size(phi)`
   pure subroutine branch_weights_with_cut(branch_weight, phi, wleb_cut, wbr, n_kept, work, keep)
      type(branch_weight_type), intent(in) :: branch_weight
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: wleb_cut
      real(wp), intent(out) :: wbr(:)
      integer, intent(out) :: n_kept
      real(wp), intent(out) :: work(:)
      logical, intent(out) :: keep(:)

      integer :: n, b, mm

      n = size(phi)
      wbr(1:n) = 1.0_wp
      n_kept = n
      if (n <= 1) return

      call branch_weight%weights(phi, wbr(1:n))
      keep(1:n) = wbr(1:n) > wleb_cut
      n_kept = count(keep(1:n))
      if (n_kept == 0) then
         n_kept = 1
         keep(1:n) = .false.
         b = maxloc(wbr(1:n), dim=1)
         keep(b) = .true.
      end if
      if (n_kept == n) return

      mm = 0
      do b = 1, n
         if (.not. keep(b)) cycle
         mm = mm + 1
         work(mm) = phi(b)
      end do
      call branch_weight%weights(work(1:n_kept), work(n_kept + 1:2*n_kept))
      mm = 0
      do b = 1, n
         if (keep(b)) then
            mm = mm + 1
            wbr(b) = work(n_kept + mm)
         else
            wbr(b) = 0.0_wp
         end if
      end do
   end subroutine branch_weights_with_cut

   !* ================================================================================= *!
   !*                         Deduplication and canonical order                         *!
   !* ================================================================================= *!

   !> Keep flags for unique feet inside the squared-distance window
   !>
   !> - a foot within `sep_cut` of an earlier kept foot is a duplicate
   !> - a foot whose `rho^2` exceeds the closest kept one by more than
   !>   `rho2_slack` carries no weight and is dropped
   !>
   !> @param[in]  pts        Feet (3, n)
   !> @param[in]  anchor     Anchor
   !> @param[in]  sep_cut    Duplicate separation (Bohr)
   !> @param[in]  rho2_slack Squared-distance window (Bohr^2)
   !> @param[out] keep       Keep flags (n)
   pure subroutine dedupe_candidates(pts, anchor, sep_cut, rho2_slack, keep)
      real(wp), intent(in) :: pts(:, :)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: sep_cut
      real(wp), intent(in) :: rho2_slack
      logical, intent(out) :: keep(:)

      integer :: i, j, m, n_kept
      real(wp) :: rho2_min, d2
      ! Kept feet so far, in run order
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

      rho2_min = huge(1.0_wp)
      do i = 1, size(pts, 2)
         if (keep(i)) rho2_min = min(rho2_min, sum((pts(:, i) - anchor)**2))
      end do
      do i = 1, size(pts, 2)
         if (.not. keep(i)) cycle
         d2 = sum((pts(:, i) - anchor)**2)
         if (d2 - rho2_min > rho2_slack) keep(i) = .false.
      end do
   end subroutine dedupe_candidates

   !> Geometry-derived branch order, so the branch index names the same surface
   !> point on every build
   !>
   !> - lexicographic on x, y, z; components closer than `tie_tol` fall through
   !>   to the next one instead of being decided by rounding noise
   !>
   !> @param[in]  pts   Feet (3, n)
   !> @param[out] order Permutation into canonical order
   pure subroutine canonical_branch_order(pts, order)
      real(wp), intent(in) :: pts(:, :)
      integer, intent(out) :: order(:)

      integer :: i, j, moved, n

      n = size(pts, 2)
      do i = 1, n
         order(i) = i
      end do
      do i = 2, n
         moved = order(i)
         j = i - 1
         do while (j >= 1)
            if (.not. precedes(moved, order(j))) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = moved
      end do

   contains

      !> Whether foot `a` sorts before foot `b`
      pure logical function precedes(a, b) result(before)
         integer, intent(in) :: a, b
         integer :: k

         before = .false.
         do k = 1, 3
            if (abs(pts(k, a) - pts(k, b)) <= tie_tol) cycle
            before = pts(k, a) < pts(k, b)
            return
         end do
      end function precedes

   end subroutine canonical_branch_order

   !* ================================================================================= *!
   !*                           Closest-point Jacobian scaling                          *!
   !* ================================================================================= *!

   !> Area scaling of the anchor-to-surface map and its switching weight
   !>
   !>     B = Q^T (w I - lambda Hess S) Q,  tau_k = Q^T t_k,  y_k = w Q B^-1 tau_k,  J = |y1 x y2|
   !>
   !> - `w_f = f_crit(||grad S||) * f_foc(beta_min(B))`; a weight that falls
   !>   below `wleb_cut` leaves `J = 1`, as the old loop did
   !>
   !> @param[in]  w            Objective weight
   !> @param[in]  g            Level set gradient at the foot
   !> @param[in]  H            Level set Hessian at the foot, packed (6)
   !> @param[in]  lambda       Multiplier
   !> @param[in]  anchor       Anchor
   !> @param[in]  owner_center Owning atom center
   !> @param[in]  f_crit       Gradient-norm switch
   !> @param[in]  f_foc        Focal switch
   !> @param[in]  wleb_cut     Weight cutoff
   !> @param[out] J            Area scaling
   !> @param[out] w_f          Switching weight
   !> @param[out] error        Singular tangent block after switching
   subroutine cp_jacobian_weight(w, g, H, lambda, anchor, owner_center, f_crit, f_foc, wleb_cut, &
                                 J, w_f, error)
      real(wp), intent(in) :: w, g(3), H(:), lambda, anchor(3), owner_center(3)
      class(moist_math_smoothstep_type), intent(in) :: f_crit
      class(moist_math_smoothstep_type), intent(in) :: f_foc
      real(wp), intent(in) :: wleb_cut
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
      if (w_f < wleb_cut) return

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

      w_f = w_f*f_foc%f0(beta2)
      if (w_f < wleb_cut) return

      if (abs(det_B) <= det_B_guard) then
         call fatal_error(error, "Tangent Jacobian matrix B is singular after switching")
         return
      end if
      Binv11 = B22/det_B
      Binv12 = -B12/det_B
      Binv22 = B11/det_B

      n_sph = anchor - owner_center
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
      call grow_r(self%rho)
      call grow_r(self%phi)
      call grow_r(self%wleb)
      call grow_r(self%wbranch)
      call grow_r(self%cpjac)
      call grow_r(self%w_f)
      call grow_r(self%r_iI)
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
      if (allocated(self%rho)) deallocate (self%rho)
      if (allocated(self%phi)) deallocate (self%phi)
      if (allocated(self%wleb)) deallocate (self%wleb)
      if (allocated(self%wbranch)) deallocate (self%wbranch)
      if (allocated(self%cpjac)) deallocate (self%cpjac)
      if (allocated(self%w_f)) deallocate (self%w_f)
      if (allocated(self%r_iI)) deallocate (self%r_iI)
      self%n = 0
      self%capacity = 0
   end subroutine rows_destroy

end module moist_cavity_drop_projection_finalize
