!> Reverse-mode (z-vector) nuclear gradient for the DROP cavity
!>
!> Contracts an already-accumulated surface adjoint directly into `dE/dR_A`,
!> at a per-point cost independent of the system size and with O(1) storage
!> (no persistent cavity derivative arrays)
!>
!> Legacy forward path in forward.f90 builds the full forward Jacobian
!> of every surface observable: `3 * n_active` seeds per grid point
!> and `O(N_sph * N_grid)` memory
!>
!>
!> A nuclear displacement reaches the surface through three routes:
!>
!>  1. Field: Nuclear dependence of the level set (field);
!>     The per-point map is linear in its seed, seeding the 13 components of the
!>     level-set jet gives adjoint weights `w_lsf0/1/2` that contract directly
!>     with the LSF's own nuclear partials
!>
!>  2. Anchor: Anchor-owner sphere-dependence (`d(anch)/dR_A = delta_{A,own}`)
!>     Adds three seeds via the objective's mixed derivative of
!>     `-d^2 phi/dr dR = alpha*I` and only affects the owner atom
!>
!>  3. Switching: The iSwiG `f_i` is function of the nuclear coordinates (and radii);
!>     It is contracted with a scalar weight
!>
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_nuclear
!$ use omp_lib, only: omp_get_thread_num
   use moist_cavity_drop_gaussian, only: iswig_workspace_type
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_surface_weights_type, seed_weight_tol
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type, seed_normal_channel, &
      & seed_jet_basis, seed_anchor, degenerate_point_error
   implicit none(type, external)

contains

   !> Contract a surface adjoint into the nuclear gradient
   !>
   !> Accumulates `dE/dR_A` for the energy whose surface adjoints `acc` holds.
   !> The result is *added* to `gradient`, so several cavities or several
   !> passes can share one accumulator.
   !>
   !> This is the reverse-mode counterpart of running `compute_gradient_drop`
   !> and contracting the resulting `*_rA` arrays; the two agree to round-off,
   !> which is what `test_cavity_drop_nuclear_adjoint` asserts.
   !>
   !> @param[in]    self     DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc      Accumulated surface-observable adjoints
   !> @param[inout] gradient Nuclear-gradient accumulator (3, nsph)
   !> @param[out]   error    Error object, allocated on failure
   module subroutine get_surface_gradient_drop(self, acc, gradient, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> Per-thread gradient buffers, summed deterministically after the region
      real(wp), allocatable :: grad_threads(:, :, :)
      !> Thread bookkeeping
      integer :: thread_slot, ithread
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort

      !> Shared per-grid point sensitivity kernel state and its response
      type(drop_seed_state_type) :: state
      !> Whether the shared prologue cleared the point for this traversal
      logical :: point_ok

      !> Grid, atom and active-slot indices
      integer :: igrid, iatom, i, n_active
      integer, allocatable :: active_idx(:)

      !> Anchor and owner sphere
      real(wp) :: anchor(3)
      integer :: owner_idx
      !> Level-set jet at the projected point
      real(wp) :: lsf1_r(3), lsf2_rr(3, 3)
      real(wp), allocatable :: lsf3_rrr(:, :, :)
      !> Jet-contracted nuclear partials of the level set, one column per active atom
      real(wp), allocatable :: vjp_pt(:, :)
      !> Objective jet at the projected point
      real(wp) :: phi1_r(3)
      !> Lagrange multiplier of the projection
      real(wp) :: lambda_val

      !> Bordered KKT system with four field seeds and three anchor seeds
      real(wp) :: kkt_rhs(4, 7)
      !> Factorization reused by every solve at this grid point
      type(drop_kkt_factor_type) :: kkt_fac

      !> Point-local level-set adjoint weights built from the 13 field seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Effective position adjoint seen by every seed
      real(wp) :: w_xyz_local(3)
      !> iSwig neighbour cache and the sparse switching rows it feeds
      type(iswig_workspace_type) :: iswig_work
      real(wp), allocatable :: swi_rows(:, :)
      real(wp) :: swi_owner_row(3), swi_f0, swi_dxi
      integer :: jj, knb

      !> Effective primitive surface adjoints
      type(drop_surface_weights_type) :: eff
      !> Timer handle
      integer :: h_sgrad

      call check_surface_adjoint(self, acc, "get_surface_gradient_drop", error)
      if (allocated(error)) return
      if (any(shape(gradient) /= [3, self%nsph])) then
         call fatal_error(error, "get_surface_gradient_drop: gradient shape mismatch")
         return
      end if
      if (self%ngrid <= 0) return

      h_sgrad = self%ctx%timer%resolve("Surface gradient", self%ctx%timer%current(), &
                                       cat_gradient)
      call self%ctx%timer%start(h_sgrad)

      !* -------------------------- Effective surface weights ------------------------- *!
      call prepare_surface_weights(self, acc, .true., eff)

      !* -------------------------------- Thread setup -------------------------------- *!
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)
      allocate (grad_threads(3, self%nsph, slots%nthreads), source=0.0_wp)

      call abort%reset()

      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& iatom, i, n_active, active_idx, state, point_ok, &
      !$omp& anchor, owner_idx, lsf1_r, lsf2_rr, lsf3_rrr, &
      !$omp& vjp_pt, phi1_r, lambda_val, &
      !$omp& kkt_rhs, kkt_fac, &
      !$omp& w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, &
      !$omp& w_xyz_local, iswig_work, swi_rows, swi_owner_row, swi_f0, swi_dxi, &
      !$omp& jj, knb)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1

      allocate (lsf3_rrr(3, 3, 3), source=0.0_wp)
      allocate (vjp_pt(3, self%nsph), source=0.0_wp)
      allocate (active_idx(self%nsph))
      call iswig_work%init(self%iswig)
      ! Sized to `nsph` rather than to the workspace capacity: `n_nb` is bounded
      ! by the atom count on either traversal, so this stays valid even if the
      ! workspace has to grow itself.
      allocate (swi_rows(3, self%nsph))

      !$omp do schedule(static, 8)
      do igrid = 1, self%ngrid
         !* -------------------------- Shared point prologue -------------------------- *!
         ! Point, jets, seed state and the solved jet and anchor seeds. Every
         ! failure path -- a latch already set, a refusing level set, a
         ! degenerate state, a singular bordered system -- has recorded itself.
         call drop_point_prologue(self, slots, thread_slot, igrid, eff%have_wk, abort, &
                                  anchor, owner_idx, lambda_val, lsf1_r, lsf2_rr, &
                                  lsf3_rrr, phi1_r, state, kkt_fac, point_ok, &
                                  kkt_rhs=kkt_rhs)
         if (.not. point_ok) cycle

         ! Outward-normal channel: the direct grad-S term rides on w_lsf1 and
         ! is picked up by the field contraction; its point-motion coupling
         ! augments the effective position weight used by every seed.
         w_lsf0_pt = 0.0_wp
         w_lsf1_pt = 0.0_wp
         w_lsf2_pt = 0.0_wp
         call seed_normal_channel(state, eff, igrid, lsf2_rr, w_lsf1_pt, w_xyz_local)

         !* -------------------- Field seeds -> level-set adjoints -------------------- *!
         call seed_jet_basis(state, eff, igrid, phi1_r, kkt_rhs, w_xyz_local, &
                             w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)

         !* -------------- Field channel: contract with nuclear partials -------------- *!
         ! The level set contracts the jet indices itself: `vjp_f1_rA` returns the
         ! nuclear-gradient row already weighted by (w_lsf0, w_lsf1, w_lsf2), so
         ! the (3, 3, 3, n_active) mixed third derivative the weights used to be
         ! folded against is never materialized -- neither here nor in the kernel.
         n_active = slots%lsf(thread_slot)%lsf%active_count()
         do i = 1, n_active
            active_idx(i) = slots%lsf(thread_slot)%lsf%active_atom(i)
         end do
         call slots%lsf(thread_slot)%lsf%vjp_f1_rA(w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, vjp_pt)
         do i = 1, n_active
            iatom = active_idx(i)
            grad_threads(:, iatom, thread_slot) = &
               grad_threads(:, iatom, thread_slot) + vjp_pt(:, i)
         end do

         !* ------------------------- Anchor channel (owner) -------------------------- *!
         call seed_anchor(state, eff, igrid, phi1_r, kkt_rhs, w_xyz_local, &
                          grad_threads(:, owner_idx, thread_slot))

         !* ------------------------- iSwig switching channel ------------------------- *!
         ! f_i depends on the nuclear geometry alone; only the owner atom and neighbours are nonzero
         if (abs(eff%w_f(igrid)) > seed_weight_tol) then
            call self%iswig%swi_collect(anchor, owner_idx, self%anchor_xi0(igrid), &
                                        swi_f0, iswig_work)
            call self%iswig%swi1_rA_sparse(iswig_work, swi_rows, swi_owner_row, swi_dxi)
            do jj = 1, iswig_work%n_nb
               knb = iswig_work%idx(jj)
               grad_threads(:, knb, thread_slot) = grad_threads(:, knb, thread_slot) &
                                                   + eff%w_f(igrid)*swi_rows(:, jj)
            end do
            grad_threads(:, owner_idx, thread_slot) = &
               grad_threads(:, owner_idx, thread_slot) + eff%w_f(igrid)*swi_owner_row
         end if

      end do
      !$omp end do

      deallocate (lsf3_rrr, vjp_pt, active_idx, swi_rows)
      call iswig_work%destroy()
      !$omp end parallel

      if (abort%requested) then
         ! An LSF failure or a KKT failure arrives as a ready-made error; a
         ! kernel degeneracy arrives as a status code that needs this routine's
         ! name to become a diagnostic.
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call degenerate_point_error("get_surface_gradient_drop", abort%status, abort%igrid, error)
         end if
         call self%ctx%timer%stop(h_sgrad)
         return
      end if

      ! Deterministic reduction: fixed thread order, independent of scheduling
      do ithread = 1, slots%nthreads
         gradient = gradient + grad_threads(:, :, ithread)
      end do

      call self%ctx%timer%stop(h_sgrad)

   end subroutine get_surface_gradient_drop

end submodule moist_cavity_drop_derivatives_nuclear
