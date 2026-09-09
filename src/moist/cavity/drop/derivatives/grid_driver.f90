!> Per-grid-point prologue shared by the DROP derivative traversals
!>
!> Three routines walk the projected grid inside an OpenMP region and open every
!> point the same way:
!>
!> * [[get_surface_gradient_drop]]  (`nuclear.f90`)
!> * [[drop_hessian_traverse]]      (`hessian_traverse.f90`), which carries both
!>   halves of the surface Hessian as two channels of one traversal
!> * [[get_surface_tangent_drop]]   (`tangent_forward.f90`)
!>
!> The common opening is: honour the abort latch, read the point, its anchor,
!> its owner sphere and its multiplier, `prepare` the thread's level set there,
!> take the level-set and objective jets, copy the grid-level scalars into the
!> kernel's seed state, derive the state, and factor the bordered KKT matrix
!> `[[H_lag, grad S], [grad S^T, 0]]` -- which is the *same* matrix in all three,
!> because all three differentiate the same stationarity condition.
!>
!> Everything past the factorization is caller specific and stays in the
!> callers. The two seams the prologue does parameterise are:
!>
!>   * `want_curvature`, forwarded to [[fill_seed_state]]. A traversal that
!>     contracts folded weights asks for whatever they carry (`eff%have_wk`);
!>     the forward tangent, and the Hessian traversal when only its
!>     adjoint-response channel is on, know their curvature channel is
!>     identically zero and ask for `.false.`
!>   * `pt%kkt_rhs`, the standard seed batch -- four level-set jet seeds and
!>     three anchor seeds. Allocated by the traversals that push the fixed
!>     16-seed basis; unallocated for [[get_surface_tangent_drop]], whose
!>     right-hand sides are one column per nuclear direction and cannot be
!>     built before the level set's directional tangents are known. That caller
!>     keeps its own solve and shares the factorization only
!>
!> `context` parameterises nothing about the computation: it is the name of the
!> API entry point the caller was reached through, and it is what turns a
!> singular bordered system at some grid point into a diagnostic the user can
!> trace back to the call they made. Every traversal above passes its own name.
!>
!> Everything the prologue produces lands in the caller's
!> [[drop_point_scratch_type]] rather than in an output argument list of its
!> own. That is what lets a traversal name one variable in its `private(...)`
!> clause instead of a dozen, and it is where the buffers that must outlive a
!> single point -- `lsf3_rrr`, `lsf4_rrrr`, the iSwiG neighbour cache -- are
!> allocated once per thread, outside the grid loop.
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_grid_driver
   use moist_cavity_drop_derivatives_kernel, only: build_seed_state, seed_state_ok
   use moist_cavity_drop_derivatives_seeds, only: seed_standard_rhs
   implicit none(type, external)

contains

   !> Prepare one grid point for a DROP derivative traversal
   !>
   !> Returns `ok = .false.` when the point must be skipped -- the latch already
   !> carries a failure, the level set refused the point, the seed state is
   !> degenerate, or the bordered system is singular. Every such exit has
   !> already recorded itself in `abort`, so the caller only has to `cycle`.
   !>
   !> `pt` is both the scratch the prologue works in and the only thing it hands
   !> back: it carries exactly the quantities the callers still read afterwards.
   !> The point itself, the level-set value, the objective value and the
   !> objective Hessian are consumed here and deliberately not stored, because
   !> nothing past the factorization asks for them.
   !>
   !> Two members of `pt` steer the prologue by their allocation status, both
   !> settled once per thread at [[point_scratch_init]]: an allocated
   !> `lsf4_rrrr` asks for the fourth level-set derivative, and an allocated
   !> `kkt_rhs` asks for the standard seed batch to be solved.
   !>
   !> `pt` is `intent(inout)` and must stay so: `intent(out)` on a derived type
   !> with allocatable components releases them on entry, which would free this
   !> thread's buffers at every grid point.
   !>
   !> @param[in]    self           DROP cavity instance
   !> @param[inout] slots          Per-thread evaluator clones (shared)
   !> @param[in]    thread_slot    Slot of the calling thread
   !> @param[in]    igrid          Grid point to prepare
   !> @param[in]    want_curvature Whether the curvature invariants are needed
   !> @param[in]    context        Calling routine, used to prefix the diagnostics
   !> @param[inout] abort          Shared failure latch of the parallel region
   !> @param[inout] pt             Point scratch of the calling thread
   !> @param[out]   ok             Whether the point may be processed further
   module subroutine drop_point_prologue(self, slots, thread_slot, igrid, &
                                         want_curvature, context, abort, pt, ok)
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
      !> Point scratch of the calling thread
      type(drop_point_scratch_type), intent(inout) :: pt
      !> Whether the point may be processed further
      logical, intent(out) :: ok

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
      pt%anchor = self%anchorxyz(:, igrid)
      pt%owner_idx = self%owner(igrid)
      pt%lambda_val = self%lambda0(igrid)

      call slots%lsf(thread_slot)%lsf%prepare(point, worker_error)

      ! The failure cannot be returned from inside the caller's worksharing
      ! construct, so park it for the post-region promotion and let the flag
      ! drain the loop. The LSF's cached derivatives are substitutes; stop
      ! before reading them.
      if (allocated(worker_error)) then
         call abort%latch_error(worker_error, igrid)
         return
      end if
      ! Every nuclear accessor of the point below shares these
      call slots%lsf(thread_slot)%lsf%cache_point_tensors()

      call slots%lsf(thread_slot)%lsf%f3_rrr(lsf0, pt%lsf1_r, pt%lsf2_rr, pt%lsf3_rrr)
      if (allocated(pt%lsf4_rrrr)) call slots%lsf(thread_slot)%lsf%f4_rrrr(pt%lsf4_rrrr)
      call slots%phi(thread_slot)%f012_r(point, pt%anchor, pt%owner_idx, phi0, &
                                         pt%phi1_r, phi2_rr)

      pt%state%lsf1_r = pt%lsf1_r
      pt%state%lsf2_rr = pt%lsf2_rr
      pt%state%lsf3_rrr = pt%lsf3_rrr
      pt%state%lambda_val = pt%lambda_val
      call fill_seed_state(self, igrid, want_curvature, pt%state)

      call build_seed_state(pt%state, self%f_crit, self%f_foc, self%f_wleb, &
                            self%param%wleb_prune_level > 0, status)
      if (status /= seed_state_ok) then
         call abort%latch_status(status, igrid)
         return
      end if

      !* ------------------------ Bordered KKT sensitivities -------------------------- *!
      ! The matrix is direction free and seed free, so every traversal factors
      ! the same thing once per grid point.
      call pt%kkt_fac%factor(phi2_rr - pt%lambda_val*pt%lsf2_rr, pt%lsf1_r, context, &
                             worker_error, igrid)
      if (allocated(worker_error)) then
         call abort%latch_error(worker_error, igrid)
         return
      end if

      if (allocated(pt%kkt_rhs)) then
         ! Which column carries which seed is not decided here: `seeds.f90` owns
         ! the layout through [[seed_rhs_column]], emits the right-hand sides
         ! that go with it, and every consumer of the solved batch --
         ! [[seed_jet_basis_apply]], [[seed_anchor_apply]], [[fill_seed_basis]]
         ! -- reads it back through the same map.
         call seed_standard_rhs(pt%lambda_val, self%param%phi_alpha, pt%kkt_rhs)
         call pt%kkt_fac%solve(pt%kkt_rhs, context, worker_error, igrid)
         if (allocated(worker_error)) then
            call abort%latch_error(worker_error, igrid)
            return
         end if
      end if

      ok = .true.

   end subroutine drop_point_prologue

end submodule moist_cavity_drop_derivatives_grid_driver
