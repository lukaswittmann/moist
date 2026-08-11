!> Reverse-mode (z-vector) nuclear gradient for the DROP cavity
!>
!> Where gradient.f90 builds the full forward Jacobian of every surface
!> observable -- `3 * n_active` seeds per grid point and `O(N_sph * N_grid)`
!> storage -- this path contracts an already-accumulated surface adjoint
!> directly into `dE/dR_A`, at a per-point cost independent of the system size
!> and with no persistent derivative arrays at all.
!>
!> A nuclear displacement reaches the surface through exactly three routes:
!>
!>  1. **Field.** The level set depends explicitly on the nuclei. Because the
!>     per-point map is linear in its seed, seeding the 13 components of the
!>     level-set jet gives adjoint weights `w_lsf0/1/2` that contract directly
!>     with the LSF's own nuclear partials.
!>
!>  2. **Anchor.** The anchor rides rigidly on its owner sphere, so
!>     `d(anchor)/dR_A = delta_{A,owner}`. This adds three seeds whose only
!>     driver is the objective's mixed derivative `-d^2 phi/dr dR = alpha*I`,
!>     and whose result lands entirely on the owner atom.
!>
!>  3. **Switching.** The iSwig overlap `f_i` is a pure function of the nuclear
!>     geometry; its gradient is already reverse-shaped, so it is contracted
!>     with a scalar weight.
!>
!> Everything else the surface depends on -- radii, the anchor Lebedev weight,
!> the sphere tangent frame -- is invariant under rigid owner motion.
submodule(moist_cavity_drop) moist_cavity_drop_nuclear_gradient
   use omp_lib, only: omp_get_thread_num
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_math_lapack_gesv, only: lapack_gesv
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_cavity_drop_adjoint_kernel, only: drop_seed_state_type, &
      & drop_seed_result_type, build_seed_state, apply_seed, compute_branch_phi_adj, &
      & seed_state_ok, seed_weight_tol, seed_status_message
   implicit none (type, external)

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
      type(lsf_thread_slot), allocatable :: lsf_threads(:)
      type(moist_cavity_drop_objective_phi_type), allocatable :: phi_threads(:)
      !> Per-thread gradient buffers, summed deterministically after the region
      real(wp), allocatable :: grad_threads(:, :, :)
      !> Thread bookkeeping
      integer :: nthreads, thread_slot, ithread
      logical :: abort_requested
      integer :: abort_status, abort_grid
      !> Per-thread LSF evaluation failure, and the one promoted for the loop
      type(error_type), allocatable :: lsf_error, abort_lsf_error

      !> Shared per-grid-point sensitivity kernel state and its response
      type(drop_seed_state_type) :: state
      type(drop_seed_result_type) :: res
      !> Degeneracy status
      integer :: status

      !> Grid, seed, atom and Cartesian indices
      integer :: igrid, ibasis, iaxis, jaxis, iatom, i, n_active
      integer, allocatable :: active_idx(:)

      !> Projected point, anchor and owner sphere
      real(wp) :: point(3), anchor(3)
      integer :: owner_idx
      !> Level-set jet at the projected point
      real(wp) :: lsf0, lsf1_r(3), lsf2_rr(3, 3)
      real(wp), allocatable :: lsf3_rrr(:, :, :)
      !> Explicit nuclear partials of the level-set jet
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :), lsf3_rr_rA(:, :, :, :)
      !> Objective jet at the projected point
      real(wp) :: phi0, phi1_r(3), phi2_rr(3, 3)
      !> Lagrange multiplier of the projection
      real(wp) :: lambda_val

      !> Bordered KKT system with four field seeds and three anchor seeds
      real(wp) :: kkt_mat_base(4, 4), kkt_mat(4, 4), kkt_rhs(4, 7)
      integer(lapack_ik) :: kkt_ipiv(4), kkt_info

      !> Seed perturbation of the level-set jet and the induced point motion
      real(wp) :: dlsf1_r(3), dlsf2_rr(3, 3)
      real(wp) :: dr_dp(3), dlambda_dp
      !> Point-local level-set adjoint weights built from the 13 field seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Accumulated adjoint contribution of one seed
      real(wp) :: contribution, g_val
      !> Outward-normal channel scratch
      real(wp) :: normal_grad(3), w_xyz_local(3), nwn
      !> iSwig switching-gradient scratch
      real(wp), allocatable :: f1_rA_pt(:, :), anchor_xi_zero(:, :)

      !> Effective primitive surface adjoints
      real(wp), allocatable :: w_xi(:), w_f_eff(:)
      real(wp), allocatable :: branch_phi_adj(:)
      logical :: have_wn, have_wk
      !> Softmax temperature
      real(wp) :: sigma_phi
      !> Timer handle
      integer :: h_sgrad

      if (.not. acc%is_initialized()) then
         call fatal_error(error, "get_surface_gradient_drop: accumulator is not initialized")
         return
      end if
      if (size(acc%w_xi) /= self%ngrid) then
         call fatal_error(error, "get_surface_gradient_drop: accumulator grid size mismatch")
         return
      end if
      if (any(shape(gradient) /= [3, self%nsph])) then
         call fatal_error(error, "get_surface_gradient_drop: gradient shape mismatch")
         return
      end if
      if (.not. allocated(self%xi0) .or. .not. allocated(self%a) .or. &
          .not. allocated(self%wleb) .or. .not. allocated(self%f)) then
         call fatal_error(error, "get_surface_gradient_drop: cavity surface data are incomplete")
         return
      end if
      if (any(abs(self%xi0) <= seed_weight_tol .and. &
              (abs(acc%w_a) > seed_weight_tol .or. abs(acc%w_w) > seed_weight_tol))) then
         call fatal_error(error, "get_surface_gradient_drop: singular derived-weight conversion")
         return
      end if
      if (self%ngrid <= 0) return

      h_sgrad = self%ctx%timer%resolve("Surface gradient", self%ctx%timer%current(), &
                                       cat_gradient)
      call self%ctx%timer%start(h_sgrad)

      !* ------------------------- Effective surface weights ------------------------- *!
      ! The area and integration-weight channels are derived:
      !   a_i = R_I^2 * f_i * wleb_i,   w_i = wleb_i,   xi_i = swx/(R_I*sqrt(wleb_i))
      ! so a_i = c*f_i/xi_i^2 and w_i = c/xi_i^2. Both fold into the width
      ! channel through d/dxi, and the area channel *additionally* folds into
      ! the switching channel through da/df = R_I^2 * wleb_i. The electronic
      ! path can skip that second fold because df/dp vanishes there; a nuclear
      ! displacement moves f, so it must not be skipped here.
      allocate (w_xi, source=acc%w_xi)
      allocate (w_f_eff, source=acc%w_f)
      where (abs(acc%w_a) > seed_weight_tol)
         w_xi = w_xi - 2.0_wp*self%a*acc%w_a/self%xi0
      end where
      where (abs(acc%w_w) > seed_weight_tol)
         w_xi = w_xi - 2.0_wp*self%wleb*acc%w_w/self%xi0
      end where
      do igrid = 1, self%ngrid
         if (abs(acc%w_a(igrid)) > seed_weight_tol) then
            w_f_eff(igrid) = w_f_eff(igrid) &
                             + acc%w_a(igrid)*self%radii(self%owner(igrid))**2*self%wleb(igrid)
         end if
      end do

      have_wn = any(abs(acc%w_n) > seed_weight_tol)
      have_wk = any(abs(acc%w_k1) > seed_weight_tol) .or. any(abs(acc%w_k2) > seed_weight_tol)

      allocate (branch_phi_adj(self%ngrid), source=0.0_wp)
      if (allocated(self%branch_count)) then
         sigma_phi = self%branch_weight%s
         call compute_branch_phi_adj(self%branch_count(1:self%ngrid), &
                                     self%anchor_id(1:self%ngrid), &
                                     self%wbranch(1:self%ngrid), &
                                     self%wleb(1:self%ngrid), &
                                     self%xi0(1:self%ngrid), &
                                     sigma_phi, w_xi, branch_phi_adj)
      end if

      !* ------------------------------- Thread setup ------------------------------- *!
      nthreads = self%ctx%get_num_threads()
      allocate (lsf_threads(nthreads))
      allocate (phi_threads(nthreads))
      do thread_slot = 1, nthreads
         allocate (lsf_threads(thread_slot)%lsf, source=self%lsf_model)
         ! The field contraction needs f3_rr_rA_screened, so the SSD cache has
         ! to be sized for third derivatives before the first %prepare
         call lsf_threads(thread_slot)%lsf%set_max_deriv(3)
         call phi_threads(thread_slot)%set_parameters(self%param)
         call phi_threads(thread_slot)%set_input(self%mol, self%radii)
      end do
      allocate (grad_threads(3, self%nsph, nthreads), source=0.0_wp)

      abort_requested = .false.
      abort_status = seed_state_ok
      abort_grid = 0

      !$omp parallel num_threads(nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& ibasis, iaxis, jaxis, iatom, i, n_active, active_idx, state, res, status, &
      !$omp& point, anchor, owner_idx, lsf0, lsf1_r, lsf2_rr, lsf3_rrr, &
      !$omp& lsf1_rA, lsf2_r_rA, lsf3_rr_rA, phi0, phi1_r, phi2_rr, lambda_val, &
      !$omp& kkt_mat_base, kkt_mat, kkt_rhs, kkt_ipiv, kkt_info, &
      !$omp& dlsf1_r, dlsf2_rr, dr_dp, dlambda_dp, &
      !$omp& w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, contribution, g_val, &
      !$omp& normal_grad, w_xyz_local, nwn, f1_rA_pt, anchor_xi_zero, lsf_error)
      thread_slot = omp_get_thread_num() + 1

      allocate (lsf3_rrr(3, 3, 3), source=0.0_wp)
      allocate (lsf1_rA(3, self%nsph), source=0.0_wp)
      allocate (lsf2_r_rA(3, 3, self%nsph), source=0.0_wp)
      allocate (active_idx(self%nsph))
      allocate (f1_rA_pt(3, self%nsph), source=0.0_wp)
      allocate (anchor_xi_zero(3, self%nsph), source=0.0_wp)

      !$omp do schedule(dynamic)
      do igrid = 1, self%ngrid
         if (abort_requested) cycle

         point = self%xyz(:, igrid)
         anchor = self%anchorxyz(:, igrid)
         owner_idx = self%owner(igrid)
         lambda_val = self%lambda0(igrid)

         call lsf_threads(thread_slot)%lsf%prepare(point, lsf_error)

         ! The failure cannot be returned from inside this worksharing construct,
         ! so park it for the post-region promotion and let the flag drain the
         ! loop. The LSF's cached derivatives are substitutes; stop before
         ! reading them.
         if (allocated(lsf_error)) then
            !$omp critical (surface_gradient_abort)
            if (.not. abort_requested) then
               abort_requested = .true.
               call move_alloc(lsf_error, abort_lsf_error)
               abort_grid = igrid
            end if
            !$omp end critical (surface_gradient_abort)
            cycle
         end if

         call lsf_threads(thread_slot)%lsf%f3_rrr_screened(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
         call lsf_threads(thread_slot)%lsf%f3_rr_rA_screened(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
         call phi_threads(thread_slot)%f012_r(point, anchor, owner_idx, phi0, phi1_r, phi2_rr)

         kkt_mat_base = 0.0_wp
         kkt_mat_base(1:3, 1:3) = phi2_rr - lambda_val*lsf2_rr
         kkt_mat_base(1:3, 4) = -lsf1_r
         kkt_mat_base(4, 1:3) = lsf1_r

         state%lsf1_r = lsf1_r
         state%lsf2_rr = lsf2_rr
         state%lsf3_rrr = lsf3_rrr
         state%lambda_val = lambda_val
         state%alpha_coeff = self%param%phi_alpha
         state%anchor = anchor
         state%owner_xyz = self%mol%xyz(:, owner_idx)
         state%anchor_wleb0 = self%anchor_wleb0(igrid)
         state%cpjac_scal0 = self%cpjac_scal0(igrid)
         state%w_f0 = self%w_f0(igrid)
         state%wbranch = self%wbranch(igrid)
         state%wleb = self%wleb(igrid)
         state%xi0 = self%xi0(igrid)
         state%want_curvature = have_wk

         call build_seed_state(state, self%f_crit, self%f_foc, self%f_wleb, &
                               self%param%wleb_prune_level > 0, status)
         if (status /= seed_state_ok) then
            !$omp critical (surface_gradient_abort)
            if (.not. abort_requested) then
               abort_requested = .true.
               abort_status = status
               abort_grid = igrid
            end if
            !$omp end critical (surface_gradient_abort)
            cycle
         end if

         ! Outward-normal channel: the direct grad-S term rides on w_lsf1 and
         ! is picked up by the field contraction; its point-motion coupling
         ! augments the effective position weight used by every seed.
         w_lsf0_pt = 0.0_wp
         w_lsf1_pt = 0.0_wp
         w_lsf2_pt = 0.0_wp
         if (have_wn) then
            nwn = dot_product(state%n_surf, acc%w_n(:, igrid))
            normal_grad = (acc%w_n(:, igrid) - state%n_surf*nwn)/state%g_norm
            w_lsf1_pt = w_lsf1_pt + normal_grad
            w_xyz_local = acc%w_xyz(:, igrid) + matmul(lsf2_rr, normal_grad)
         else
            w_xyz_local = acc%w_xyz(:, igrid)
         end if

         !* ------------------------- Bordered KKT sensitivities ------------------------- *!
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
         kkt_mat = kkt_mat_base
         call lapack_gesv(4_lapack_ik, 7_lapack_ik, kkt_mat, 4_lapack_ik, &
                          kkt_ipiv, kkt_rhs, 4_lapack_ik, kkt_info)
         if (kkt_info /= 0_lapack_ik) then
            !$omp critical (surface_gradient_abort)
            if (.not. abort_requested) then
               abort_requested = .true.
               abort_status = -1
               abort_grid = igrid
            end if
            !$omp end critical (surface_gradient_abort)
            cycle
         end if

         !* ---------------------- Field seeds -> level-set adjoints ---------------------- *!
         do ibasis = 1, 13
            dlsf1_r = 0.0_wp
            dlsf2_rr = 0.0_wp
            if (ibasis <= 4) then
               if (ibasis > 1) dlsf1_r(ibasis - 1) = 1.0_wp
               dr_dp = kkt_rhs(1:3, ibasis)
               dlambda_dp = kkt_rhs(4, ibasis)
            else
               iaxis = (ibasis - 5)/3 + 1
               jaxis = mod(ibasis - 5, 3) + 1
               dlsf2_rr(iaxis, jaxis) = 1.0_wp
               dr_dp = 0.0_wp
               dlambda_dp = 0.0_wp
            end if

            call apply_seed(state, dlsf1_r, dlsf2_rr, dr_dp, dlambda_dp, res)

            contribution = dot_product(w_xyz_local, dr_dp) + w_xi(igrid)*res%dxi
            if (abs(branch_phi_adj(igrid)) > seed_weight_tol) then
               contribution = contribution + branch_phi_adj(igrid)*dot_product(phi1_r, dr_dp)
            end if
            if (have_wk) then
               contribution = contribution + acc%w_k1(igrid)*res%dk1 + acc%w_k2(igrid)*res%dk2
            end if

            if (ibasis == 1) then
               w_lsf0_pt = w_lsf0_pt + contribution
            else if (ibasis <= 4) then
               w_lsf1_pt(ibasis - 1) = w_lsf1_pt(ibasis - 1) + contribution
            else
               w_lsf2_pt(iaxis, jaxis) = w_lsf2_pt(iaxis, jaxis) + contribution
            end if
         end do

         !* ------------------ Field channel: contract with nuclear partials ------------------ *!
         n_active = lsf_threads(thread_slot)%lsf%active_count()
         do i = 1, n_active
            active_idx(i) = lsf_threads(thread_slot)%lsf%active_atom(i)
         end do
         do i = 1, n_active
            iatom = active_idx(i)
            do iaxis = 1, 3
               g_val = w_lsf0_pt*lsf1_rA(iaxis, iatom) &
                       + dot_product(w_lsf1_pt, lsf2_r_rA(:, iaxis, iatom)) &
                       + sum(w_lsf2_pt*lsf3_rr_rA(:, :, iaxis, iatom))
               grad_threads(iaxis, iatom, thread_slot) = &
                  grad_threads(iaxis, iatom, thread_slot) + g_val
            end do
         end do

         !* ---------------------------- Anchor channel (owner) ---------------------------- *!
         do iaxis = 1, 3
            dr_dp = kkt_rhs(1:3, 4 + iaxis)
            dlambda_dp = kkt_rhs(4, 4 + iaxis)
            dlsf1_r = 0.0_wp
            dlsf2_rr = 0.0_wp

            call apply_seed(state, dlsf1_r, dlsf2_rr, dr_dp, dlambda_dp, res)

            contribution = dot_product(w_xyz_local, dr_dp) + w_xi(igrid)*res%dxi
            if (abs(branch_phi_adj(igrid)) > seed_weight_tol) then
               ! phi = 0.5*alpha*|r - anchor|^2, so at fixed r the owner's rigid
               ! motion contributes -phi1_r on top of the point-motion term
               contribution = contribution + branch_phi_adj(igrid) &
                              *(dot_product(phi1_r, dr_dp) - phi1_r(iaxis))
            end if
            if (have_wk) then
               contribution = contribution + acc%w_k1(igrid)*res%dk1 + acc%w_k2(igrid)*res%dk2
            end if

            grad_threads(iaxis, owner_idx, thread_slot) = &
               grad_threads(iaxis, owner_idx, thread_slot) + contribution
         end do

         !* --------------------------- iSwig switching channel --------------------------- *!
         ! f_i depends on the nuclear geometry alone; swi1_rA already returns
         ! the full (3, nsph) gradient, so this channel costs the same in both
         ! modes. anchor_xi has no nuclear dependence, matching gradient.f90.
         if (abs(w_f_eff(igrid)) > seed_weight_tol) then
            f1_rA_pt = self%iswig%swi1_rA(anchor, owner_idx, self%anchor_xi0(igrid), &
                                          anchor_xi_zero)
            grad_threads(:, :, thread_slot) = grad_threads(:, :, thread_slot) &
                                              + w_f_eff(igrid)*f1_rA_pt
         end if

      end do
      !$omp end do

      deallocate (lsf3_rrr, lsf1_rA, lsf2_r_rA, active_idx, f1_rA_pt, anchor_xi_zero)
      if (allocated(lsf3_rr_rA)) deallocate (lsf3_rr_rA)
      !$omp end parallel

      if (abort_requested) then
         if (allocated(abort_lsf_error)) then
            call move_alloc(abort_lsf_error, error)
         else if (abort_status == -1) then
            call fatal_error(error, "get_surface_gradient_drop: KKT sensitivity solve failed")
         else
            call degenerate_gradient_point_error(abort_status, abort_grid, error)
         end if
         call self%ctx%timer%stop(h_sgrad)
         return
      end if

      ! Deterministic reduction: fixed thread order, independent of scheduling
      do ithread = 1, nthreads
         gradient = gradient + grad_threads(:, :, ithread)
      end do

      call self%ctx%timer%stop(h_sgrad)

   end subroutine get_surface_gradient_drop

   !> Report a degenerate projected point from the sensitivity kernel
   !>
   !> @param[in]  status  One of the `seed_state_*` codes
   !> @param[in]  igrid   Offending grid point
   !> @param[out] error   Error object
   subroutine degenerate_gradient_point_error(status, igrid, error)
      !> Degeneracy status
      integer, intent(in) :: status
      !> Offending grid point
      integer, intent(in) :: igrid
      !> Error object
      type(error_type), allocatable, intent(out) :: error

      !> Rendered grid-point index
      character(len=32) :: idx

      write (idx, '(i0)') igrid
      call fatal_error(error, "[Error] get_surface_gradient_drop: "// &
                       seed_status_message(status)//" at grid point "//trim(idx))

   end subroutine degenerate_gradient_point_error

end submodule moist_cavity_drop_nuclear_gradient
