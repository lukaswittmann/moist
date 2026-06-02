!> DROP projection and mapped Jacobian quadrature weight computation
submodule(moist_cavity_drop) moist_cavity_drop_projection
   use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_wtime
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none

contains

   !* ================================================================================= *!
   !*             Memory management utilities for projection-coupled arrays             *!
   !* ================================================================================= *!

   !> Compute next storage capacity using +10% growth steps
   !>
   !> TODO: This is not optimal yet
   module pure integer function projection_grow_capacity(current_capacity, required_capacity) result(new_capacity)
      integer, intent(in) :: current_capacity
      integer, intent(in) :: required_capacity

      integer :: increment

      new_capacity = max(1, current_capacity)
      do while (new_capacity < required_capacity)
         increment = max(1, ceiling(0.1_wp*real(new_capacity, wp)))
         new_capacity = new_capacity + increment
      end do
   end function projection_grow_capacity

   !> Ensure all projection-coupled arrays share at least the requested capacity
   !>
   !> TODO: Error propagation is missing
   module subroutine ensure_projection_capacity(self, new_capacity)
      class(cavity_type_drop), intent(inout) :: self
      integer, intent(in) :: new_capacity

      !> Grow 2D real arrays
      call grow_array(self%xyz, 3, new_capacity, fill_value=0.0_wp)
      call grow_array(self%anchorxyz, 3, new_capacity, fill_value=0.0_wp)
      call grow_array(self%normal0, 3, new_capacity, fill_value=0.0_wp)

      !> Grow 1D real arrays
      call grow_array(self%wleb, new_capacity, fill_value=0.0_wp)
      call grow_array(self%anchor_wleb0, new_capacity, fill_value=0.0_wp)
      call grow_array(self%lambda0, new_capacity, fill_value=0.0_wp)
      call grow_array(self%iswig_f0, new_capacity, fill_value=1.0_wp)
      call grow_array(self%f, new_capacity, fill_value=1.0_wp)
      call grow_array(self%anchor_xi0, new_capacity, fill_value=0.0_wp)
      call grow_array(self%rho, new_capacity, fill_value=0.0_wp)
      call grow_array(self%r_iI0, new_capacity, fill_value=0.0_wp)
      call grow_array(self%wbranch, new_capacity, fill_value=1.0_wp)
      call grow_array(self%phi0, new_capacity, fill_value=0.0_wp)

      !> Grow 1D integer arrays
      call grow_array(self%owner, new_capacity, fill_value=0)
      call grow_array(self%branch, new_capacity, fill_value=1)
      call grow_array(self%anchor_id, new_capacity, fill_value=0)
      call grow_array(self%branch_count, new_capacity, fill_value=1)

      !> Grow 1D logical arrays
      call grow_array(self%converged, new_capacity, fill_value=.false.)
   end subroutine ensure_projection_capacity

   !* ================================================================================= *!
   !*                                     Projection                                    *!
   !* ================================================================================= *!

   !> Project all grid points onto SDF surface
   !>
   !> Uses [[drop_projector_type]] to minimize the objective function subject to S=0
   !>
   !> @param[inout] self Cavity with grid to project
   module subroutine project_all_points(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      type(prettyprinter) :: pp
      type(prettylistprinter) :: plp
      integer :: i, idx, ithread
      integer :: iend, ibeg, nloc
      integer :: nmax_anchor, nout, n_branch
      integer :: n_branched_anchor, n_branched_points
      integer :: nbranch_min, nbranch_max
      integer :: nthreads, thread_slot
      integer :: number_base
      integer :: local_branched_anchor, local_branched_points
      integer :: local_nbranch_min, local_nbranch_max
      integer :: local_done, est_done, next_progress_pct
      integer :: n_append_fail, n_proj_fail
      real(wp) :: progress_last_time, now_time, progress_elapsed, progress_rate, progress_eta
      real(wp) :: wall_start, wall_end
      logical :: append_ok
      logical :: trigger_time, trigger_pct
      integer, allocatable :: thread_nout(:), thread_branched_anchor(:), thread_branched_points(:)
      integer, allocatable :: thread_nbranch_min(:), thread_nbranch_max(:)
      type(drop_projector_type), allocatable :: projectors(:)
      type(projection_buffer_type), allocatable :: proj_buffers(:)
      type(projection_workspace_type), allocatable :: works(:)

      ! nmax_anchor = input anchors; nout = output points (as points can branch)
      nmax_anchor = self%nmax
      nout = 0

      nthreads = max(1, omp_get_max_threads())

      ! Private per-thread projector, output buffer and scratch workspace
      allocate (projectors(nthreads))
      allocate (proj_buffers(nthreads))
      allocate (works(nthreads))

      ! Per-thread tallies, reduced serially below (for now)
      allocate (thread_nout(nthreads), source=0)
      allocate (thread_branched_anchor(nthreads), source=0)
      allocate (thread_branched_points(nthreads), source=0)
      allocate (thread_nbranch_min(nthreads), source=huge(1))
      allocate (thread_nbranch_max(nthreads), source=0)

      ! Buffer preallocation: expected share + 10% slack to avoid mid-loop regrow
      nloc = max(1, ceiling(1.1_wp*real(max(1, nmax_anchor/max(1, nthreads)), wp)))
      n_append_fail = 0
      n_proj_fail = 0

      !*------------------------- Project all points------------------------- *!

      if (self%verbosity >= 2) then
         plp = new_prettylistprinter([12, 8, 12, 12], &
                                     [character(len=8) :: "Point", "% done", "Elapsed (s)", "Left (s)"], &
                                     unit=output_unit)
         call plp%blank()
         call plp%header("PROJECTOR")
         call plp%blank()
         call plp%print_header()
         call plp%separator()
      end if

      wall_start = omp_get_wtime()

      !$omp parallel num_threads(nthreads) default(shared) private(thread_slot, i, n_branch, append_ok, &
      !$omp& local_branched_anchor, local_branched_points, local_nbranch_min, local_nbranch_max, &
      !$omp& local_done, est_done, next_progress_pct, progress_last_time, now_time, progress_elapsed, &
      !$omp& progress_rate, progress_eta, trigger_time, trigger_pct)
      thread_slot = omp_get_thread_num() + 1

      ! First-touch init inside the region so allocations are core-local
      call projectors(thread_slot)%destroy()

      ! init() sets solver tolerances and branching cutoffs
      call projectors(thread_slot)%init(self%param, self%lsf_model, &
                                        branch_sep_cut=self%param%branch_sep_cut, &
                                        branch_rho_cut=self%param%branch_rho_cut, &
                                        verbosity=self%verbosity, &
                                        debug=self%debug, tol=self%param%proj_tol, maxiter=self%param%proj_maxiter)
      call proj_buffers(thread_slot)%init(nloc)
      call works(thread_slot)%init()
      !init_primitives() binds the objective/LSF to this molecule, radii and screening grid
      call projectors(thread_slot)%init_primitives(self%mol, self%radii, self%mol_cell_grid)

      local_branched_anchor = 0
      local_branched_points = 0
      local_nbranch_min = huge(1)
      local_nbranch_max = 0
      local_done = 0
      progress_last_time = wall_start
      next_progress_pct = 10

      !$omp do schedule(dynamic)
      do i = 1, nmax_anchor

         ! Reset workspace size (keeps allocation) so no leftover branches survive
         call works(thread_slot)%clear()

         ! f below cutoff; downstream weight filtering will remove that point
         if (self%f(i) < self%param%wleb_cut) then
            if (self%verbosity >= 3) then
               !$omp critical (projection_warning_print)
               write (output_unit, '(x,a,i0,a)') &
                  "Skipping projection for gridpoint ", i, " (f below cutoff)"
               !$omp end critical (projection_warning_print)
            end if
            call works(thread_slot)%set_single(self%anchorxyz(:, i), .true.)
         else
            block
               type(error_type), allocatable :: proj_error

               ! Actually project point (can return multiple branches)
               call projectors(thread_slot)%project_point( &
                  anchor=self%anchorxyz(:, i), &
                  owner=self%owner(i), &
                  initial_guess=self%anchorxyz(:, i), &
                  index=i, &
                  error=proj_error, &
                  proj_level=self%param%proj_level, &
                  n_points=n_branch, &
                  work=works(thread_slot) &
                  )

               if (allocated(proj_error)) then
                  !$omp atomic update
                  n_proj_fail = n_proj_fail + 1
                  !$omp critical (projection_warning_print)
                  print '(a,i0)', "Warning: Projection failed for gridpoint ", i
                  print '(a)', proj_error%message
                  print '(a)', "Debug info:"
                  print '(a,i0)', "  gridpoint ID: ", i
                  print '(a,3(es12.4))', "  anchor: ", self%anchorxyz(:, i)
                  !$omp end critical (projection_warning_print)
                  ! On failure keep the anchor but flag it not converged so it is
                  ! later assigned f=0 and excluded from the quadrature.
                  call works(thread_slot)%set_single(self%anchorxyz(:, i), .false.)
               end if
            end block
         end if

         ! Branches actually produced for this anchor (1 normally, >1 if branched)
         n_branch = works(thread_slot)%size()
         call works(thread_slot)%reserve(n_branch)

         ! Branch weights are computed once on the final surviving set in compute_branch_weights after filter
         works(thread_slot)%branch_weights(1:n_branch) = 1.0_wp

         if (n_branch > 1) then
            local_branched_anchor = local_branched_anchor + 1
            local_branched_points = local_branched_points + n_branch
            local_nbranch_min = min(local_nbranch_min, n_branch)
            local_nbranch_max = max(local_nbranch_max, n_branch)
         end if

         ! Copy this anchor's branches plus its per-anchor scalars into the thread-local buffer
         call proj_buffers(thread_slot)%add_workspace( &
            work=works(thread_slot), &
            anchor_xyz=self%anchorxyz(:, i), &
            owner=self%owner(i), &
            wleb=self%wleb(i), &
            anchor_wleb0=self%anchor_wleb0(i), &
            iswig_f0=self%iswig_f0(i), &
            f=self%f(i), &
            anchor_xi0=self%anchor_xi0(i), &
            anchor_id=self%anchor_id(i), &
            ok=append_ok)
         if (.not. append_ok) then
            !$omp atomic update
            n_append_fail = n_append_fail + 1
         end if

         ! Progress is reported only by thread 1, extrapolating its own count to all threads (est_done)
         local_done = local_done + 1
         if (thread_slot == 1) then
            if (self%verbosity >= 2) then
               est_done = min(nmax_anchor, local_done*nthreads)
               now_time = omp_get_wtime()

               trigger_time = (now_time - progress_last_time) >= 5.0_wp
               trigger_pct = nmax_anchor > 0 .and. (100*est_done >= next_progress_pct*nmax_anchor)

               if (trigger_time .or. trigger_pct) then
                  progress_elapsed = now_time - wall_start
                  progress_rate = 0.0_wp
                  progress_eta = 0.0_wp
                  if (progress_elapsed > 0.0_wp) then
                     progress_rate = real(est_done, wp)/progress_elapsed
                     if (progress_rate > tiny(1.0_wp)) then
                        progress_eta = real(max(0, nmax_anchor - est_done), wp)/progress_rate
                     end if
                  end if

                  !$omp critical (projection_warning_print)
                  call plp%begin_row()
                  call plp%add(est_done)
                  call plp%add(100.0_wp*real(est_done, wp)/real(max(1, nmax_anchor), wp), fmt='f6.2')
                  call plp%add(progress_elapsed, fmt='f8.2')
                  call plp%add(progress_eta, fmt='f8.2')
                  call plp%end_row()
                  !$omp end critical (projection_warning_print)

                  progress_last_time = now_time
                  do while (next_progress_pct <= 100 .and. 100*est_done >= next_progress_pct*nmax_anchor)
                     next_progress_pct = next_progress_pct + 10
                  end do
               end if
            end if
         end if

      end do ! i = 1, nmax_anchor
      !$omp end do

      ! Publish this thread's branch stats into the shared arrays (one slot each)
      thread_branched_anchor(thread_slot) = local_branched_anchor
      thread_branched_points(thread_slot) = local_branched_points
      thread_nbranch_min(thread_slot) = local_nbranch_min
      thread_nbranch_max(thread_slot) = local_nbranch_max
      !$omp end parallel

      wall_end = omp_get_wtime()
      if (self%verbosity >= 2) then
         call plp%separator()
         call plp%begin_row()
         call plp%add(nmax_anchor)
         call plp%add("")
         call plp%add(wall_end - wall_start, fmt='f8.2')
         call plp%add(0.0_wp, fmt='f8.2')
         call plp%end_row()
      end if

      if (n_append_fail > 0) then
         call fatal_error(error, "Error: Invalid branch payload in projection buffer append.")
         return
      end if

      ! Total output count is now known: sum each thread buffer's size
      do ithread = 1, nthreads
         thread_nout(ithread) = proj_buffers(ithread)%size()
         nout = nout + thread_nout(ithread)
      end do

      ! Grow the cavity's coupled arrays to fit nout
      call ensure_projection_capacity(self, max(1, nout))
      self%rho = 0.0_wp
      self%r_iI0 = 0.0_wp
      self%converged = .false.
      self%normal0 = 0.0_wp
      self%lambda0 = 0.0_wp
      self%branch = 1
      self%anchor_id = 0
      self%branch_count = 1
      self%wbranch = 1.0_wp
      self%phi0 = 0.0_wp

      ! Concatenate thread buffers into the flat cavity arrays, in thread order (deterministic ordering)
      idx = 0
      if (nout > 0) then
         do ithread = 1, nthreads
            nloc = thread_nout(ithread)
            if (nloc <= 0) cycle

            ! Contiguous output slice [ibeg:iend] for this thread's points
            ibeg = idx + 1
            iend = idx + nloc

            self%xyz(:, ibeg:iend) = proj_buffers(ithread)%xyz(:, 1:nloc)
            self%anchorxyz(:, ibeg:iend) = proj_buffers(ithread)%anchorxyz(:, 1:nloc)
            self%normal0(:, ibeg:iend) = proj_buffers(ithread)%normal0(:, 1:nloc)

            self%wleb(ibeg:iend) = proj_buffers(ithread)%wleb(1:nloc)
            self%anchor_wleb0(ibeg:iend) = proj_buffers(ithread)%anchor_wleb0(1:nloc)
            self%lambda0(ibeg:iend) = proj_buffers(ithread)%lambda0(1:nloc)
            self%iswig_f0(ibeg:iend) = proj_buffers(ithread)%iswig_f0(1:nloc)
            self%f(ibeg:iend) = proj_buffers(ithread)%f(1:nloc)
            self%anchor_xi0(ibeg:iend) = proj_buffers(ithread)%anchor_xi0(1:nloc)
            self%rho(ibeg:iend) = proj_buffers(ithread)%rho(1:nloc)
            self%wbranch(ibeg:iend) = proj_buffers(ithread)%wbranch(1:nloc)
            self%phi0(ibeg:iend) = proj_buffers(ithread)%phi0(1:nloc)

            self%owner(ibeg:iend) = proj_buffers(ithread)%owner(1:nloc)
            self%branch(ibeg:iend) = proj_buffers(ithread)%branch(1:nloc)
            self%anchor_id(ibeg:iend) = proj_buffers(ithread)%anchor_id(1:nloc)
            self%branch_count(ibeg:iend) = proj_buffers(ithread)%branch_count(1:nloc)

            self%converged(ibeg:iend) = proj_buffers(ithread)%converged(1:nloc)

            idx = iend
         end do

         ! Per-point finalize: converged points get their owner-sphere distance
         do idx = 1, nout
            if (self%converged(idx)) then
               self%r_iI0(idx) = norm2(self%xyz(:, idx) - self%mol%xyz(:, self%owner(idx)))
            else
               self%r_iI0(idx) = 0.0_wp
               self%f(idx) = 0.0_wp
            end if
         end do
      end if

      self%ngrid = nout

      ! Build numbering array for all points:
      !   This is a unique global id = anchor_id + number_base*(branch-1) and is (has been) very helpful
      !   for debugging and tracking grid points through displacements
      if (allocated(self%numbering)) deallocate (self%numbering)
      allocate (self%numbering(nout), source=-1)
      number_base = max(1, self%nsph*self%param%num_leb)
      do idx = 1, nout
         self%numbering(idx) = self%anchor_id(idx) + number_base*(self%branch(idx) - 1)
      end do

      n_branched_anchor = sum(thread_branched_anchor)
      n_branched_points = sum(thread_branched_points)

      if (self%verbosity >= 1) then
         pp = new_prettyprinter(unit=output_unit)
         call pp%blank()
         call pp%push('Projection results:')
         call pp%kv('Failed projections', n_proj_fail)
         call pp%kv('Workspace append fails', n_append_fail)
         call pp%kv('Total points', nout)
         call pp%kv('Anchors with branches', n_branched_anchor)
         call pp%kv('New branched points', n_branched_points)
         call pp%pop()
      end if

   end subroutine project_all_points

   !* ================================================================================= *!
   !*                           Closest-point Jacobian scaling                          *!
   !* ================================================================================= *!

   !> Compute closest-point Jacobian scaling
   !>
   !> Computes the area scaling factor J_i for each projected grid point by
   !> evaluating the Jacobian of the projection map from the anchor sphere to
   !> the SDF surface. Works entirely in the tangent basis as
   !>
   !>   Q = [q1, q2]                     surface tangent frame  (from n = g/|g|)
   !>   B = Q^T A Q                      tangent-restricted KKT matrix  (2x2 symmetric)
   !>   tau_k = Q^T t_k                  sphere tangent vectors projected into surface tangent plane
   !>   y_k = alpha * Q * B^{-1} * tau_k
   !>   J_i = |y1 x y2|
   !>
   !> A smooth eigenvalue-based switch turns off points as the smallest tangent
   !> eigenvalue approaches zero or becomes negative (focal point)
   !>
   !> @param[inout] self Cavity with projected grid points
   module subroutine compute_cp_jacobian_scaling(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: igrid, nthreads, thread_slot
      real(wp) :: proj_point(3), anchor_point(3), lambda_val
      real(wp) :: lsf0
      real(wp), allocatable :: lsf1_r_threads(:, :), lsf2_rr_threads(:, :, :)
      real(wp) :: A(3, 3), g_vec(3), g_norm_sq, g_norm
      real(wp) :: alpha_coeff
      type(lsf_thread_slot), allocatable :: lsf_threads(:)
      ! Surface tangent frame Q = [q1, q2] from n = g/|g|
      real(wp) :: n_surf(3), q1(3), q2(3)

      ! Tangent-restricted KKT matrix B = Q^T A Q (2x2 symmetric)
      real(wp) :: B11, B12, B22
      ! Analytic 2x2 eigenvalues of B
      real(wp) :: tr_B, det_B, disc, sqrt_disc
      real(wp) :: beta1, beta2   ! eigenvalues (beta1 >= beta2)
      real(wp) :: lambda_switch_i
      real(wp), parameter :: det_B_guard = 1.0e-30_wp

      ! Unregularized inverse of B
      real(wp) :: Binv11, Binv12, Binv22

      ! Sphere tangent frame and tangent-space projection
      real(wp) :: n_sph(3), t1(3), t2(3)
      real(wp) :: tau1(2), tau2(2)  ! tau_k = Q^T t_k
      real(wp) :: w1(2), w2(2)     ! w_k = Binv_reg * tau_k

      ! Jacobian computation
      real(wp) :: y1(3), y2(3), cross_prod(3), J_i

      logical :: abort_requested

      ! Tangent-restricted KKT diagnostics (debug only)
      logical :: do_diag
      integer :: n_diag_points
      type(prettyprinter) :: pp
      type(prettylistprinter) :: plp_diag
      logical, allocatable :: diag_mask(:)
      integer :: c_critical, c_warning, c_safe

      abort_requested = .false.

      ! Initialize SSD systems and thread-local SSD evaluators
      nthreads = max(1, omp_get_max_threads())
      allocate (lsf_threads(nthreads))
      allocate (lsf1_r_threads(3, nthreads), source=0.0_wp)
      allocate (lsf2_rr_threads(3, 3, nthreads), source=0.0_wp)
      do thread_slot = 1, nthreads
         allocate (lsf_threads(thread_slot)%lsf, source=self%lsf_model)
      end do

      ! Initialize Jacobian scaling array
      if (allocated(self%cpjac_scal0)) deallocate (self%cpjac_scal0)
      allocate (self%cpjac_scal0(self%ngrid), source=1.0_wp)

      if (allocated(self%w_f0)) deallocate (self%w_f0)
      allocate (self%w_f0(self%ngrid), source=1.0_wp)

      ! Coefficient alpha = 0.5 * w_a from the Jacobian formula:
      ! J_P = alpha * A^{-1} * (I - g*g^T*A^{-1}/(g^T*A^{-1}*g))
      alpha_coeff = self%param%phi_alpha

      ! Loop over all grid points to compute Jacobian scaling
      !$omp parallel num_threads(nthreads) default(shared) private(thread_slot, igrid, proj_point, &
      !$omp& anchor_point, lambda_val, lsf0, A, g_vec, g_norm_sq, g_norm, n_surf, q1, q2, &
      !$omp& B11, B12, B22, tr_B, det_B, disc, sqrt_disc, beta1, beta2, lambda_switch_i, &
      !$omp& Binv11, Binv12, Binv22, n_sph, t1, t2, tau1, tau2, w1, w2, y1, y2, cross_prod, J_i)
      thread_slot = omp_get_thread_num() + 1

      !$omp do schedule(dynamic)
      do igrid = 1, self%ngrid
         !$omp cancellation point do
         if (abort_requested) cycle

         ! Skip if below switching cutoff.
         if (self%f(igrid) < self%param%wleb_cut) then
            self%w_f0(igrid) = 0.0_wp
            cycle
         end if

         proj_point = self%xyz(:, igrid)
         anchor_point = self%anchorxyz(:, igrid)
         lambda_val = self%lambda0(igrid)

         ! Compute SSD on-the-fly for this point
         call lsf_threads(thread_slot)%lsf%prepare(proj_point)

         ! Compute lsf derivatives
         call lsf_threads(thread_slot)%lsf%f012_r_screened(lsf0, &
                                                           lsf1_r_threads(:, thread_slot), lsf2_rr_threads(:, :, thread_slot))

         g_vec = lsf1_r_threads(:, thread_slot)
         g_norm_sq = dot_product(g_vec, g_vec)
         g_norm = sqrt(g_norm_sq)

         ! Compute ||g||-based switching function
         self%w_f0(igrid) = self%f_crit%f0(g_norm)

         ! Skip if weight below cutoff
         if (self%w_f0(igrid) < self%param%wleb_cut) cycle

         ! Build A = (w_a/2) * I - lambda * Hess(SDF)
         A(1, 1) = alpha_coeff - lambda_val*lsf2_rr_threads(1, 1, thread_slot)
         A(1, 2) = -lambda_val*lsf2_rr_threads(1, 2, thread_slot)
         A(1, 3) = -lambda_val*lsf2_rr_threads(1, 3, thread_slot)
         A(2, 1) = -lambda_val*lsf2_rr_threads(2, 1, thread_slot)
         A(2, 2) = alpha_coeff - lambda_val*lsf2_rr_threads(2, 2, thread_slot)
         A(2, 3) = -lambda_val*lsf2_rr_threads(2, 3, thread_slot)
         A(3, 1) = -lambda_val*lsf2_rr_threads(3, 1, thread_slot)
         A(3, 2) = -lambda_val*lsf2_rr_threads(3, 2, thread_slot)
         A(3, 3) = alpha_coeff - lambda_val*lsf2_rr_threads(3, 3, thread_slot)

         ! Surface tangent frame Q = [q1, q2] from n = g/|g|
         n_surf = g_vec/g_norm
         call setup_tangent_frame(n_surf, q1, q2)

         ! Tangent-restricted KKT matrix B = Q^T A Q (2x2 symmetric)
         B11 = dot_product(q1, matmul(A, q1))
         B12 = dot_product(q1, matmul(A, q2))
         B22 = dot_product(q2, matmul(A, q2))

         ! Analytic 2x2 eigenvalues of B
         ! beta_{1,2} = tr(B)/2 +/- sqrt(tr(B)^2/4 - det(B))
         tr_B = B11 + B22
         det_B = B11*B22 - B12*B12
         disc = 0.25_wp*tr_B*tr_B - det_B
         disc = max(disc, 0.0_wp)
         sqrt_disc = sqrt(disc)
         beta1 = 0.5_wp*tr_B + sqrt_disc ! larger eigenvalue
         beta2 = 0.5_wp*tr_B - sqrt_disc ! smaller eigenvalue

         lambda_switch_i = beta2
         self%w_f0(igrid) = self%w_f0(igrid)*self%f_foc%f0(lambda_switch_i)
         if (self%w_f0(igrid) < self%param%wleb_cut) cycle

         if (abs(det_B) <= det_B_guard) then
            !$omp critical (compute_cpjac_abort)
            if (.not. abort_requested) then
               abort_requested = .true.
               call fatal_error(error, "[Error] Tangent Jacobian matrix B is singular after switching")
            end if
            !$omp end critical (compute_cpjac_abort)
            !$omp cancel do
            cycle
         end if

         ! Direct inverse of B
         Binv11 = B22/det_B
         Binv12 = -B12/det_B
         Binv22 = B11/det_B

         ! Sphere tangent frame (for the anchor point)
         n_sph = anchor_point - self%mol%xyz(:, self%owner(igrid))
         call setup_tangent_frame(n_sph, t1, t2)

         ! Project sphere tangent vectors into surface tangent plane
         ! tau_k = Q^T t_k  (2-vectors)
         tau1(1) = dot_product(q1, t1)
         tau1(2) = dot_product(q2, t1)
         tau2(1) = dot_product(q1, t2)
         tau2(2) = dot_product(q2, t2)

         ! Compute w_k = Binv * tau_k
         w1(1) = Binv11*tau1(1) + Binv12*tau1(2)
         w1(2) = Binv12*tau1(1) + Binv22*tau1(2)
         w2(1) = Binv11*tau2(1) + Binv12*tau2(2)
         w2(2) = Binv12*tau2(1) + Binv22*tau2(2)

         ! Lift back to 3D: y_k = alpha * Q * w_k
         y1 = alpha_coeff*(w1(1)*q1 + w1(2)*q2)
         y2 = alpha_coeff*(w2(1)*q1 + w2(2)*q2)

         ! Cross product: y1 x y2
         cross_prod(1) = y1(2)*y2(3) - y1(3)*y2(2)
         cross_prod(2) = y1(3)*y2(1) - y1(1)*y2(3)
         cross_prod(3) = y1(1)*y2(2) - y1(2)*y2(1)

         ! Area scaling J_i = |y1 x y2|
         J_i = sqrt(cross_prod(1)**2 + cross_prod(2)**2 + cross_prod(3)**2)
         self%cpjac_scal0(igrid) = J_i

      end do ! igrid
      !$omp end do
      !$omp end parallel

      if (allocated(error)) return
      if (abort_requested) then
         call fatal_error(error, "Jacobian scaling computation aborted. (unreachable in normal execution)")
         return
      end if

   end subroutine compute_cp_jacobian_scaling

   !* ================================================================================= *!
   !*                              Gaussian surface charges                             *!
   !* ================================================================================= *!

   !> Recompute gaussian surface charge widths with projected wleb (after CP Jacobian scaling)
   !> This is the "true" xi used for cpcm energy and gradients.
   module subroutine compute_gaussians(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: i

      if (allocated(self%xi0)) deallocate (self%xi0)
      allocate (self%xi0(self%ngrid), source=0.0_wp)

      do i = 1, self%ngrid
         self%xi0(i) = self%iswig%xi0(self%owner(i), self%wleb(i))
      end do

   end subroutine compute_gaussians

end submodule moist_cavity_drop_projection
