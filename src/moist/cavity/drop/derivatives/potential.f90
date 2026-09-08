!> Reverse-mode surface -> level set adjoint contractions for the DROP cavity.
!>
!> Provides the variational cavity response (Fock) infrastucture
!>
!> These routines map per-point surface adjoint weights (Gaussian width,
!> integration weight, area, switch factor, projected position, and normal) onto
!> adjoint weights of the level set function value/gradient/Hessian
!>
!> The per-grid point sensitivity kernel is shared with the nuclear path
!> in [[moist_cavity_drop_derivatives_kernel]]
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_potential
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type, &
      & drop_point_scratch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type
   use moist_cavity_drop_derivatives_seeds, only: seed_normal_channel, seed_jet_basis
   implicit none(type, external)

contains

   !> Map accumulated surface adjoints into the generic response container
   !>
   !> @param[inout] self     DROP cavity instance
   !> @param[in]    acc      Accumulated surface-observable adjoints
   !> @param[inout] response Response accumulator receiving the LSF channels
   !> @param[out]   error    Error object
   module subroutine get_surface_response_drop(self, acc, response, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(inout) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Response accumulator receiving the LSF channels
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)

      allocate (w0(self%ngrid), w1(3, self%ngrid), w2(3, 3, self%ngrid))
      call self%contract_surface_lsf_weights(acc, w0, w1, w2, error)
      if (allocated(error)) return

      if (.not. allocated(response%lsf%w_value)) then
         allocate (response%lsf%w_value(self%ngrid), source=0.0_wp)
         allocate (response%lsf%w_gradient(3, self%ngrid), source=0.0_wp)
         allocate (response%lsf%w_hessian(3, 3, self%ngrid), source=0.0_wp)
      else if (size(response%lsf%w_value) /= self%ngrid) then
         call fatal_error(error, "DROP surface response grid-size mismatch")
         return
      end if
      response%lsf%w_value = response%lsf%w_value + w0
      response%lsf%w_gradient = response%lsf%w_gradient + w1
      response%lsf%w_hessian = response%lsf%w_hessian + w2

   end subroutine get_surface_response_drop

   !> Contract per-grid surface adjoint weights to LSF value/gradient/Hessian adjoints
   !>
   !> Implements the DROP reverse-mode chain rule: a perturbation $p$
   !> of the level set function (value $S$, gradient $\nabla S$, Hessian $\nabla^2 S$) moves
   !> the projected surface point and every weight derived from it.
   !>
   !> This routine rewrites the upstream surface adjoint as the LSF-local adjoint,
   !> enforcing for every perturbation $p$ the identity
   !>
   !> $$
   !> \sum_i \Big[ w^{\xi}_i \, \frac{\partial \xi_i}{\partial p}
   !>            + \mathbf{w}^{\mathrm{xyz}}_i \cdot \frac{\partial \mathbf{r}_i}{\partial p} \Big]
   !> = \sum_i \Big[ w^{S}_i \, \frac{\partial S_i}{\partial p}
   !>            + \mathbf{w}^{S_r}_i \cdot \frac{\partial (\nabla S_i)}{\partial p}
   !>            + \sum_{a,b} w^{S_{rr}}_{ab,i} \, \frac{\partial (\nabla^2 S_i)_{ab}}{\partial p} \Big],
   !> $$
   !>
   !> where $w^{\xi}$ = `acc%w_xi`, $\mathbf{w}^{\mathrm{xyz}}$ = `acc%w_xyz`, and
   !> $w^{S}, \mathbf{w}^{S_r}, w^{S_{rr}}$ = `w_lsf0`, `w_lsf1`, `w_lsf2`. The
   !> outward-normal (`acc%w_n`) and principal-curvature (`acc%w_k1`, `acc%w_k2`)
   !> channels enter the same left-hand sum through their own $\partial/\partial p$
   !> sensitivities and are folded into the same `w_lsf` weights. The derived
   !> area (`acc%w_a`) and integration-weight (`acc%w_w`) channels are folded into
   !> the Gaussian-width channel via $a_i = c f_i/\xi_i^2$ and $w_i = c/\xi_i^2$.
   !>
   !> The switching factor $f_i$ is an anchor-only iSwig overlap, so $\partial
   !> f_i/\partial p = 0$ for a level-set perturbation and `acc%w_f` does not
   !> contribute here. It does contribute to the nuclear gradient, which is why
   !> the area channel is folded into the width channel only.
   !>
   !> The per-point opening is [[drop_point_prologue]]'s, exactly as on the four
   !> parallel traversals, so this electronic path cannot drift away from them.
   !> It runs the prologue on a one-slot [[drop_worker_slots_type]] and a local
   !> latch, and because the loop is serial it reads that latch itself rather
   !> than draining the grid: an LSF or KKT failure is returned immediately, and
   !> a degenerate seed state clears the latch and skips its point alone.
   !>
   !> The prologue solves the full seven-column seed batch where this routine
   !> needs only the four jet columns. The three extra right-hand sides are
   !> three more columns of one `getrs` on an already factored 4x4, which is not
   !> worth a second batch layout to avoid; [[seed_jet_basis]] reads columns 1-4
   !> and ignores the rest.
   !>
   !> @param[in]  self    DROP cavity instance (must hold a projected grid)
   !> @param[in]  acc     Accumulated surface-observable adjoints
   !> @param[out] w_lsf0  Adjoint weights for LSF values S_i (ngrid)
   !> @param[out] w_lsf1  Adjoint weights for LSF gradients S_r_i (3, ngrid)
   !> @param[out] w_lsf2  Adjoint weights for LSF Hessians S_rr_i (3, 3, ngrid)
   !> @param[out] error   Error object, allocated on failure (KKT sensitivity solve)
   module subroutine contract_surface_lsf_weights(self, acc, w_lsf0, w_lsf1, w_lsf2, error)
      !> DROP cavity instance (must hold a projected grid)
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Adjoint weights for LSF values S_i (ngrid)
      real(wp), intent(out) :: w_lsf0(:)
      !> Adjoint weights for LSF gradients S_r_i (3, ngrid)
      real(wp), intent(out) :: w_lsf1(:, :)
      !> Adjoint weights for LSF Hessians S_rr_i (3, 3, ngrid)
      real(wp), intent(out) :: w_lsf2(:, :, :)
      !> Error object, allocated on failure (KKT sensitivity solve)
      type(error_type), allocatable, intent(out) :: error

      !> Level-set clone and objective used to rebuild the per-point jet; one
      !> slot, because this traversal is serial
      type(drop_worker_slots_type) :: slots
      !> Per-thread state of the grid point being opened
      type(drop_point_scratch_type) :: pt
      !> Failure latch of the prologue; read and cleared by this loop itself
      type(drop_abort_latch_type) :: abort
      !> Whether the shared prologue cleared the point
      logical :: point_ok
      !> Grid index
      integer :: igrid
      !> Point-local level-set adjoints built from the 13 jet seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Folded surface adjoints and the branch objective adjoint
      type(drop_surface_weights_type) :: eff
      real(wp) :: w_xyz_local(3)

      call check_surface_adjoint(self, acc, "contract_surface_lsf_weights", error)
      if (allocated(error)) return

      ! fold_switching = .false.: the electronic degrees of freedom leave the
      ! switching factor f untouched, so the area channel's da/df term is
      ! identically zero here. The nuclear path passes .true.
      call prepare_surface_weights(self, acc, .false., eff)

      ! One slot rather than the context's team: nothing below is parallel, and
      ! a clone carries the level set's screened-derivative cache with it.
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii, &
                      nthreads=1)
      call pt%init(self%nsph, self%iswig, want_seed_batch=.true.)
      call abort%reset()

      w_lsf0 = 0.0_wp
      w_lsf1 = 0.0_wp
      w_lsf2 = 0.0_wp

      do igrid = 1, self%ngrid
         !* -------------------------- Shared point prologue -------------------------- *!
         ! Point, jets, seed state and the solved jet and anchor seeds -- the
         ! same opening the four parallel traversals get.
         call drop_point_prologue(self, slots, 1, igrid, eff%have_wk, &
                                  "contract_surface_lsf_weights", abort, pt, point_ok)
         if (.not. point_ok) then
            ! The latch exists for a worksharing construct that cannot return.
            ! This loop can, so it reads the latch itself: a level-set refusal
            ! or a singular bordered system arrives as an error object and is
            ! returned at once, while a degenerate seed state carries only a
            ! status code and skips its own point. That one has to clear the
            ! latch, or the next prologue call would see a failure still
            ! requested and drain the rest of the grid.
            if (allocated(abort%error)) then
               call move_alloc(abort%error, error)
               return
            end if
            call abort%reset()
            cycle
         end if

         ! Fold an optional outward-normal adjoint weight into the field channels:
         ! the direct grad-S contribution normal_grad = P_tan(w_n)/|grad S| enters
         ! w_lsf1 at the fixed projected point, and its point-motion coupling
         ! H @ normal_grad augments the effective position weight below.
         !
         ! This write happens only once the point is known to be usable, so a
         ! rejected point never leaves a half-contracted weight behind.
         w_lsf0_pt = 0.0_wp
         w_lsf1_pt = 0.0_wp
         w_lsf2_pt = 0.0_wp
         call seed_normal_channel(pt%state, eff, igrid, pt%lsf2_rr, w_lsf1_pt, w_xyz_local)

         ! The 13 jet seeds share the point's one factorization. Only the value
         ! (column 1) and the three gradient directions (columns 2-4) move the
         ! point; the nine Hessian perturbations have rhs = 0 and hence
         ! dr/dp = 0, dlambda/dp = 0.
         call seed_jet_basis(pt%state, eff, igrid, pt%phi1_r, pt%kkt_rhs, w_xyz_local, &
                             w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)

         w_lsf0(igrid) = w_lsf0(igrid) + w_lsf0_pt
         w_lsf1(:, igrid) = w_lsf1(:, igrid) + w_lsf1_pt
         w_lsf2(:, :, igrid) = w_lsf2(:, :, igrid) + w_lsf2_pt
      end do

      call pt%destroy()

   end subroutine contract_surface_lsf_weights

end submodule moist_cavity_drop_derivatives_potential
