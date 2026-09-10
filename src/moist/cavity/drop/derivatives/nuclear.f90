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
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type, &
      & drop_point_scratch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type, seed_weight_tol
   use moist_cavity_drop_derivatives_seeds, only: seed_normal_channel, &
      & seed_jet_basis, seed_anchor
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

      !> Per-thread state of the grid point being opened: the jets, the seed
      !> state, the bordered factorization, the solved seed batch and the
      !> buffers they are read into
      type(drop_point_scratch_type) :: pt
      !> Whether the shared prologue cleared the point for this traversal
      logical :: point_ok

      !> Grid, atom and active-slot indices
      integer :: igrid, iatom, i

      !> Point-local level-set adjoint weights built from the 13 field seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Effective position adjoint seen by every seed
      real(wp) :: w_xyz_local(3)
      !> Sparse switching rows of the owner sphere
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
      !$omp& pt, point_ok, iatom, i, &
      !$omp& w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, &
      !$omp& w_xyz_local, swi_owner_row, swi_f0, swi_dxi, &
      !$omp& jj, knb)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1

      call pt%init(self%nsph, self%iswig, want_seed_batch=.true., want_vjp=.true.)

      !$omp do schedule(static, 8)
      do igrid = 1, self%ngrid
         !* -------------------------- Shared point prologue -------------------------- *!
         ! Point, jets, seed state and the solved jet and anchor seeds. Every
         ! failure path -- a latch already set, a refusing level set, a
         ! degenerate state, a singular bordered system -- has recorded itself.
         call drop_point_prologue(self, slots, thread_slot, igrid, eff%have_wk, &
                                  "get_surface_gradient_drop", abort, pt, point_ok)
         if (.not. point_ok) cycle

         ! Outward-normal channel: the direct grad-S term rides on w_lsf1 and
         ! is picked up by the field contraction; its point-motion coupling
         ! augments the effective position weight used by every seed.
         w_lsf0_pt = 0.0_wp
         w_lsf1_pt = 0.0_wp
         w_lsf2_pt = 0.0_wp
         call seed_normal_channel(pt%state, eff, igrid, pt%lsf2_rr, w_lsf1_pt, w_xyz_local)

         !* -------------------- Field seeds -> level-set adjoints -------------------- *!
         call seed_jet_basis(pt%state, eff, igrid, pt%phi1_r, pt%kkt_rhs, w_xyz_local, &
                             w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)

         !* -------------- Field channel: contract with nuclear partials -------------- *!
         ! The level set contracts the jet indices itself: `vjp_f1_rA` returns the
         ! nuclear-gradient row already weighted by (w_lsf0, w_lsf1, w_lsf2), so
         ! the (3, 3, 3, n_active) mixed third derivative the weights used to be
         ! folded against is never materialized -- neither here nor in the kernel.
         pt%n_active = slots%lsf(thread_slot)%lsf%active_count()
         do i = 1, pt%n_active
            pt%active_idx(i) = slots%lsf(thread_slot)%lsf%active_atom(i)
         end do
         call slots%lsf(thread_slot)%lsf%vjp_f1_rA(w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, pt%vjp_pt)
         do i = 1, pt%n_active
            iatom = pt%active_idx(i)
            grad_threads(:, iatom, thread_slot) = &
               grad_threads(:, iatom, thread_slot) + pt%vjp_pt(:, i)
         end do

         !* ------------------------- Anchor channel (owner) -------------------------- *!
         call seed_anchor(pt%state, eff, igrid, pt%phi1_r, pt%kkt_rhs, w_xyz_local, &
                          grad_threads(:, pt%owner_idx, thread_slot))

         !* ------------------------- iSwig switching channel ------------------------- *!
         ! f_i depends on the nuclear geometry alone; only the owner atom and neighbours are nonzero
         if (abs(eff%w_f(igrid)) > seed_weight_tol) then
            call self%iswig%swi_collect(pt%anchor, pt%owner_idx, self%anchor_xi0(igrid), &
                                        swi_f0, pt%iswig_work)
            call self%iswig%swi1_rA_sparse(pt%iswig_work, pt%swi_rows, swi_owner_row, swi_dxi)
            do jj = 1, pt%iswig_work%n_nb
               knb = pt%iswig_work%idx(jj)
               grad_threads(:, knb, thread_slot) = grad_threads(:, knb, thread_slot) &
                                                   + eff%w_f(igrid)*pt%swi_rows(:, jj)
            end do
            grad_threads(:, pt%owner_idx, thread_slot) = &
               grad_threads(:, pt%owner_idx, thread_slot) + eff%w_f(igrid)*swi_owner_row
         end if

      end do
      !$omp end do

      call pt%destroy()
      !$omp end parallel

      if (abort%requested) then
         call abort%raise("get_surface_gradient_drop", error)
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
