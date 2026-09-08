!> Per-grid-point prologue shared by the DROP derivative traversals
!>
!> Four routines walk the projected grid inside an OpenMP region and open every
!> point the same way:
!>
!> * [[get_surface_gradient_drop]]         (`nuclear.f90`)
!> * [[get_surface_hessian_fixed_drop]]    (`hessian_fixed.f90`)
!> * [[get_surface_hessian_response_drop]] (`hessian_response.f90`)
!> * [[get_surface_tangent_drop]]          (`tangent_forward.f90`)
!>
!> The common opening is: honour the abort latch, read the point, its anchor,
!> its owner sphere and its multiplier, `prepare` the thread's level set there,
!> take the level-set and objective jets, copy the grid-level scalars into the
!> kernel's seed state, derive the state, and factor the bordered KKT matrix
!> `[[H_lag, grad S], [grad S^T, 0]]` -- which is the *same* matrix in all four,
!> because all four differentiate the same stationarity condition.
!>
!> Everything past the factorization is caller specific and stays in the
!> callers. The two seams the prologue does parameterise are:
!>
!>   * `want_curvature`, forwarded to [[fill_seed_state]]. The two reverse
!>     traversals ask for whatever their folded weights carry (`eff%have_wk`);
!>     the response and forward-tangent halves know their curvature channel is
!>     identically zero and ask for `.false.`
!>   * `kkt_rhs`, the standard seed batch -- four level-set jet seeds and three
!>     anchor seeds. Present for the three traversals that push the fixed
!>     16-seed basis; absent for [[get_surface_tangent_drop]], whose right-hand
!>     sides are one column per nuclear direction and cannot be built before
!>     the level set's directional tangents are known. That caller keeps its
!>     own solve and shares the factorization only
!>
!> `context` parameterises nothing about the computation: it is the name of the
!> API entry point the caller was reached through, and it is what turns a
!> singular bordered system at some grid point into a diagnostic the user can
!> trace back to the call they made. Every traversal above passes its own name.
!>
!> `lsf3_rrr` and `lsf4_rrrr` stay caller-owned buffers rather than becoming
!> locals here: they are allocated once per thread outside the grid loop, and
!> an automatic array would move that allocation into the hot loop.
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_grid_driver
   use moist_cavity_drop_derivatives_kernel, only: build_seed_state, seed_state_ok
   implicit none(type, external)

contains

   !> Prepare one grid point for a DROP derivative traversal
   !>
   !> Returns `ok = .false.` when the point must be skipped -- the latch already
   !> carries a failure, the level set refused the point, the seed state is
   !> degenerate, or the bordered system is singular. Every such exit has
   !> already recorded itself in `abort`, so the caller only has to `cycle`.
   !>
   !> The outputs are exactly the quantities the callers still read after the
   !> prologue. The point itself, the level-set value, the objective value and
   !> the objective Hessian are consumed here and deliberately not handed back,
   !> which is what lets the callers drop them from their `private(...)` lists.
   !>
   !> @param[in]    self           DROP cavity instance
   !> @param[inout] slots          Per-thread evaluator clones (shared)
   !> @param[in]    thread_slot    Slot of the calling thread
   !> @param[in]    igrid          Grid point to prepare
   !> @param[in]    want_curvature Whether the curvature invariants are needed
   !> @param[in]    context        Calling routine, used to prefix the diagnostics
   !> @param[inout] abort          Shared failure latch of the parallel region
   !> @param[out]   anchor         Anchor of the grid point
   !> @param[out]   owner_idx      Owner sphere of the anchor
   !> @param[out]   lambda_val     Lagrange multiplier of the projection
   !> @param[out]   lsf1_r         Level-set gradient at the projected point
   !> @param[out]   lsf2_rr        Level-set Hessian at the projected point
   !> @param[inout] lsf3_rrr       Caller-owned third-derivative buffer
   !> @param[out]   phi1_r         Objective gradient at the projected point
   !> @param[inout] state          Seed state
   !> @param[out]   kkt_fac        Factorized bordered KKT system
   !> @param[out]   ok             Whether the point may be processed further
   !> @param[inout] lsf4_rrrr      Caller-owned fourth-derivative buffer
   !> @param[out]   kkt_rhs        Solved standard jet and anchor seeds
   module subroutine drop_point_prologue(self, slots, thread_slot, igrid, &
                                         want_curvature, context, abort, anchor, &
                                         owner_idx, lambda_val, lsf1_r, lsf2_rr, &
                                         lsf3_rrr, phi1_r, state, kkt_fac, ok, &
                                         lsf4_rrrr, kkt_rhs)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type), intent(inout) :: slots
      !> Slot of the calling thread
      integer, intent(in) :: thread_slot
      !> Grid point to prepare
      integer, intent(in) :: igrid
      !> Whether the curvature invariants are needed
      logical, intent(in) :: want_curvature
      !> Calling routine, so a failure names the entry point the user called
      character(len=*), intent(in) :: context
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type), intent(inout) :: abort
      !> Anchor of the grid point
      real(wp), intent(out) :: anchor(3)
      !> Owner sphere of the anchor
      integer, intent(out) :: owner_idx
      !> Lagrange multiplier of the projection
      real(wp), intent(out) :: lambda_val
      !> Level-set gradient and Hessian at the projected point
      real(wp), intent(out) :: lsf1_r(3), lsf2_rr(3, 3)
      !> Caller-owned third-derivative buffer
      real(wp), intent(inout) :: lsf3_rrr(3, 3, 3)
      !> Objective gradient at the projected point
      real(wp), intent(out) :: phi1_r(3)
      !> Shared per-grid point sensitivity kernel state
      type(drop_seed_state_type), intent(inout) :: state
      !> Factorization reused by every solve at this grid point
      type(drop_kkt_factor_type), intent(out) :: kkt_fac
      !> Whether the point may be processed further
      logical, intent(out) :: ok
      !> Caller-owned fourth-derivative buffer
      real(wp), intent(inout), optional :: lsf4_rrrr(3, 3, 3, 3)
      !> Standard batch of four jet seeds and three anchor seeds, already solved
      real(wp), intent(out), optional :: kkt_rhs(4, 7)

      !> Projected point
      real(wp) :: point(3)
      !> Level-set and objective values, and the objective Hessian; read by the
      !> jet accessors and not needed past the factorization
      real(wp) :: lsf0, phi0, phi2_rr(3, 3)
      !> Degeneracy status
      integer :: status
      !> Failure on its way to the latch
      type(error_type), allocatable :: worker_error

      ok = .false.
      if (abort%requested) return

      point = self%xyz(:, igrid)
      anchor = self%anchorxyz(:, igrid)
      owner_idx = self%owner(igrid)
      lambda_val = self%lambda0(igrid)

      call slots%lsf(thread_slot)%lsf%prepare(point, worker_error)

      ! The failure cannot be returned from inside the caller's worksharing
      ! construct, so park it for the post-region promotion and let the flag
      ! drain the loop. The LSF's cached derivatives are substitutes; stop
      ! before reading them.
      if (allocated(worker_error)) then
         call abort%latch_error(worker_error, igrid)
         return
      end if

      call slots%lsf(thread_slot)%lsf%f3_rrr(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      if (present(lsf4_rrrr)) call slots%lsf(thread_slot)%lsf%f4_rrrr(lsf4_rrrr)
      call slots%phi(thread_slot)%f012_r(point, anchor, owner_idx, phi0, phi1_r, phi2_rr)

      state%lsf1_r = lsf1_r
      state%lsf2_rr = lsf2_rr
      state%lsf3_rrr = lsf3_rrr
      state%lambda_val = lambda_val
      call fill_seed_state(self, igrid, want_curvature, state)

      call build_seed_state(state, self%f_crit, self%f_foc, self%f_wleb, &
                            self%param%wleb_prune_level > 0, status)
      if (status /= seed_state_ok) then
         call abort%latch_status(status, igrid)
         return
      end if

      !* ------------------------ Bordered KKT sensitivities -------------------------- *!
      ! The matrix is direction free and seed free, so every traversal factors
      ! the same thing once per grid point.
      call kkt_fac%factor(phi2_rr - lambda_val*lsf2_rr, lsf1_r, context, &
                          worker_error, igrid)
      if (allocated(worker_error)) then
         call abort%latch_error(worker_error, igrid)
         return
      end if

      if (present(kkt_rhs)) then
         ! Columns 1-4 are the level-set value and gradient seeds; the nine
         ! Hessian seeds have a zero right-hand side. Columns 5-7 are the
         ! anchor seeds: moving the owner rigidly leaves the field untouched
         ! and drives the system through -d^2 phi/dr dR_owner = +alpha*I.
         kkt_rhs = 0.0_wp
         kkt_rhs(4, 1) = -1.0_wp
         kkt_rhs(1, 2) = lambda_val
         kkt_rhs(2, 3) = lambda_val
         kkt_rhs(3, 4) = lambda_val
         kkt_rhs(1, 5) = self%param%phi_alpha
         kkt_rhs(2, 6) = self%param%phi_alpha
         kkt_rhs(3, 7) = self%param%phi_alpha
         call kkt_fac%solve(kkt_rhs, context, worker_error, igrid)
         if (allocated(worker_error)) then
            call abort%latch_error(worker_error, igrid)
            return
         end if
      end if

      ok = .true.

   end subroutine drop_point_prologue

end submodule moist_cavity_drop_derivatives_grid_driver
