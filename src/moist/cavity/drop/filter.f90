!> DROP grid filtering routines.
submodule(moist_cavity_drop) moist_cavity_drop_filter
   use mctc_env_accuracy, only: wp
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none

contains

   !> Compact all per-grid arrays using a keep mask
   !>
   !> This helper centralizes the actual filtering logic so different
   !> filtering stages only need to define the keep criterion.
   !>
   !> @param[inout] self   Cavity data to compact
   !> @param[in]    nold   Number of active grid points before filtering
   !> @param[in]    keep   Logical keep mask of length `nold`
   !> @param[in]    nvalid Number of surviving grid points
   subroutine compact_grid_arrays(self, nold, keep, nvalid)
      class(cavity_type_drop), intent(inout) :: self
      integer, intent(in) :: nold, nvalid
      logical, intent(in) :: keep(:)
      integer :: i, number_base
      integer, allocatable :: expected_numbering(:)

      if (allocated(self%numbering)) then
         call filter_array(self%numbering, nold, keep, nvalid)
      end if

      call filter_array(self%xyz, nold, keep, nvalid)
      call filter_array(self%anchorxyz, nold, keep, nvalid)
      call filter_array(self%normal0, nold, keep, nvalid)

      call filter_array(self%wleb, nold, keep, nvalid)
      call filter_array(self%anchor_wleb0, nold, keep, nvalid)
      call filter_array(self%lambda0, nold, keep, nvalid)
      call filter_array(self%iswig_f0, nold, keep, nvalid)
      call filter_array(self%f, nold, keep, nvalid)
      call filter_array(self%anchor_xi0, nold, keep, nvalid)
      call filter_array(self%rho, nold, keep, nvalid)
      call filter_array(self%r_iI0, nold, keep, nvalid)
      call filter_array(self%wbranch, nold, keep, nvalid)
      call filter_array(self%phi0, nold, keep, nvalid)
      call filter_array(self%cpjac_scal0, nold, keep, nvalid)
      call filter_array(self%w_f0, nold, keep, nvalid)

      call filter_array(self%owner, nold, keep, nvalid)
      call filter_array(self%branch, nold, keep, nvalid)
      call filter_array(self%anchor_id, nold, keep, nvalid)
      call filter_array(self%branch_count, nold, keep, nvalid)

      call filter_array(self%converged, nold, keep, nvalid)

      self%ngrid = nvalid

   end subroutine compact_grid_arrays

   !> Remove grid points below switching cutoff using the current point weights
   !>
   !> @param[inout] self Cavity with projected grid
   module subroutine filter_arrays(self, name, error)
      class(cavity_type_drop), intent(inout) :: self
      character(len=*), intent(in) :: name
      type(error_type), allocatable, intent(out) :: error
      type(prettyprinter) :: pp

      integer :: ncur, nvalid, nremoved
      logical, allocatable :: keep(:)

      ncur = self%ngrid
      allocate (keep(ncur), source=.false.)

      !> Remove based on the current accumulated quadrature weights.
      keep = self%wleb*self%f > self%param%wleb_cut

      nvalid = count(keep)
      nremoved = ncur - nvalid
      call compact_grid_arrays(self, ncur, keep, nvalid)

      if (self%verbosity >= 2) then
         pp = new_prettyprinter(unit=output_unit)
         call pp%blank()
         call pp%kv(name//' removed', nremoved, 'points')
      end if

   end subroutine filter_arrays

   !> Compute branch weights on the final surviving set of grid points
   !> and fold them into wleb.
   !>
   !> For each contiguous anchor group (group_size > 1) this routine:
   !>   1. Computes softmax over the full set of surviving phi values
   !>   2. Marks any sibling whose weight is below `wleb_cut` for
   !>      removal (sets its wleb to zero so the next filter_arrays
   !>      call drops it)
   !>   3. Re-computes softmax restricted to the kept siblings, so
   !>      their weights once again sum to 1.
   !>   4. Writes wbranch, multiplies wleb by wbranch, and sets
   !>      branch_count to the kept-sibling count
   !>
   !> Singleton groups (group_size == 1) keep their placeholder
   !> wbranch = 1.0 from projection and are skipped
   !>
   !> The caller is responsible for a follow-up filter_arrays call to
   !> compact the grid after any siblings have been marked for removal
   module subroutine compute_branch_weights(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      integer :: igroup_start, igroup_end, group_size, m, im_grid
      integer :: owner_idx, n_kept, mm
      real(wp) :: anchor_wleb0_val
      real(wp), allocatable :: phi_local(:), weights_new(:)
      real(wp), allocatable :: phi_kept(:), weights_kept(:)
      logical, allocatable :: keep(:)
      type(prettylistprinter) :: plp_branch, plp_dist
      type(prettyprinter) :: pp
      integer :: mi, mj
      real(wp) :: d_ij
      integer, allocatable :: dist_widths(:)
      character(len=16), allocatable :: dist_headers(:)
      character(len=16) :: row_label

      ! Pre-scan / post-loop statistics
      integer :: igrp_s, igrp_e, grp_sz
      integer :: pre_n_branched_groups, pre_n_branched_points
      integer :: pre_nbranch_max
      real(wp) :: pre_avg_branches
      integer :: post_n_multi_groups, post_n_kept_points, post_n_dropped_points

      ! Kept points only across anchors that *remain* multi-branch (n_kept > 1);
      ! post_n_kept_points includes singletons collapsed from originally
      ! branched groups, so it is not the right denominator/numerator for
      ! "branched points (kept)" or "branches/anchor (avg)"
      integer :: post_n_multi_kept
      ! Anchors that started branched but collapsed to a single sibling.
      integer :: post_n_collapsed_groups
      integer :: post_nbranch_max
      real(wp) :: post_avg_kept

      ! Pre-scan: walk the grid in anchor-id groups and accumulate
      ! per-branched-group statistics without touching any state
      pre_n_branched_groups = 0
      pre_n_branched_points = 0
      pre_nbranch_max = 0
      igrp_s = 1
      do while (igrp_s <= self%ngrid)
         if (self%branch_count(igrp_s) <= 1) then
            igrp_s = igrp_s + 1
            cycle
         end if
         igrp_e = igrp_s
         do while (igrp_e < self%ngrid)
            if (self%anchor_id(igrp_e + 1) /= self%anchor_id(igrp_s)) exit
            igrp_e = igrp_e + 1
         end do
         grp_sz = igrp_e - igrp_s + 1
         pre_n_branched_groups = pre_n_branched_groups + 1
         pre_n_branched_points = pre_n_branched_points + grp_sz
         pre_nbranch_max = max(pre_nbranch_max, grp_sz)
         igrp_s = igrp_e + 1
      end do
      if (pre_n_branched_groups > 0) then
         pre_avg_branches = real(pre_n_branched_points, wp)/real(pre_n_branched_groups, wp)
      else
         pre_avg_branches = 0.0_wp
      end if

      if (self%verbosity >= 2) then
         pp = new_prettyprinter(unit=output_unit)
         call pp%blank()
         call pp%push('Branch weighting (initial):')
         call pp%kv('Branched anchors', pre_n_branched_groups)
         call pp%kv('Branched points', pre_n_branched_points)
         call pp%kv('Branches/anchor (avg)', pre_avg_branches)
         call pp%kv('Branches/anchor (max)', pre_nbranch_max)
         call pp%pop()
      end if

      ! Post-loop accumulators populated while processing each group.
      post_n_multi_groups = 0
      post_n_multi_kept = 0
      post_n_collapsed_groups = 0
      post_n_kept_points = 0
      post_n_dropped_points = 0
      post_nbranch_max = 0

      ! Per-branched-group detail table
      if (self%verbosity >= 3 .and. pre_n_branched_groups > 0) then
         plp_branch = new_prettylistprinter([12, 10, 10, 16, 16, 16], &
                                            [character(len=16) :: 'anchor_id', 'owner', 'branch', 'phi', 'distance', 'weight'], &
                                            unit=output_unit, fmt_len=16)
      end if

      igroup_start = 1
      do while (igroup_start <= self%ngrid)
         if (self%branch_count(igroup_start) <= 1) then
            igroup_start = igroup_start + 1
            cycle
         end if

         ! Extend the group while anchor_id stays the same
         igroup_end = igroup_start
         do while (igroup_end < self%ngrid)
            if (self%anchor_id(igroup_end + 1) /= self%anchor_id(igroup_start)) exit
            igroup_end = igroup_end + 1
         end do
         group_size = igroup_end - igroup_start + 1

         ! Initial softmax over the whole group to decide which siblings are meaningful
         owner_idx = self%owner(igroup_start)
         anchor_wleb0_val = self%anchor_wleb0(igroup_start)
         allocate (phi_local(group_size), weights_new(group_size), keep(group_size))
         do m = 1, group_size
            phi_local(m) = self%phi0(igroup_start + m - 1)
         end do
         call self%branch_weight%weights(phi_local, weights_new)

         ! Keep siblings whose weight would stay above wleb_cut
         keep = weights_new > self%param%wleb_cut
         n_kept = count(keep)

         if (n_kept == 0) then
            ! Safety fallback: keep the single strongest sibling if the
            ! threshold would have thrown the entire group away
            n_kept = 1
            keep = .false.
            keep(maxloc(weights_new, dim=1)) = .true.
         end if

         if (n_kept < group_size) then
            ! Re-softmax over kept phi values so the retained branch weights sum to 1 exactly
            allocate (phi_kept(n_kept), weights_kept(n_kept))
            mm = 0
            do m = 1, group_size
               if (keep(m)) then
                  mm = mm + 1
                  phi_kept(mm) = phi_local(m)
               end if
            end do
            call self%branch_weight%weights(phi_kept, weights_kept)
            mm = 0
            do m = 1, group_size
               if (keep(m)) then
                  mm = mm + 1
                  weights_new(m) = weights_kept(mm)
               else
                  weights_new(m) = 0.0_wp
               end if
            end do
            deallocate (phi_kept, weights_kept)
         end if

         do m = 1, group_size
            im_grid = igroup_start + m - 1
            if (keep(m)) then
               self%wbranch(im_grid) = weights_new(m)
               self%wleb(im_grid) = self%wleb(im_grid)*weights_new(m)
               self%branch_count(im_grid) = n_kept
            else
               ! Mark for removal by the follow-up filter_arrays call
               self%wbranch(im_grid) = 0.0_wp
               self%wleb(im_grid) = 0.0_wp
               self%branch_count(im_grid) = n_kept
            end if
         end do

         ! Post-loop accumulation: tally kept and dropped survivors per group
         post_n_kept_points = post_n_kept_points + n_kept
         post_n_dropped_points = post_n_dropped_points + (group_size - n_kept)
         post_nbranch_max = max(post_nbranch_max, n_kept)
         if (n_kept > 1) then
            post_n_multi_groups = post_n_multi_groups + 1
            post_n_multi_kept = post_n_multi_kept + n_kept
         else
            post_n_collapsed_groups = post_n_collapsed_groups + 1
         end if

         ! Per-group diagnostics
         if (self%verbosity >= 3) then
            call plp_branch%blank()
            call plp_branch%print_header()
            call plp_branch%separator()
            do m = 1, group_size
               call plp_branch%begin_row()
               call plp_branch%add(self%anchor_id(igroup_start))
               call plp_branch%add(owner_idx)
               call plp_branch%add(m)
               call plp_branch%add(phi_local(m), fmt='es12.4')
               call plp_branch%add(self%rho(igroup_start + m - 1), fmt='f14.8')
               call plp_branch%add(weights_new(m), fmt='es12.4')
               call plp_branch%end_row()
            end do

            ! Pairwise distance matrix between projected sibling points
            allocate (dist_widths(group_size + 1))
            allocate (dist_headers(group_size + 1))
            dist_widths(1) = 10
            dist_widths(2:) = 16
            dist_headers(1) = 'i \ j'
            do mj = 1, group_size
               write (dist_headers(mj + 1), '(i0)') mj
            end do
            plp_dist = new_prettylistprinter(dist_widths, dist_headers, &
                                             unit=output_unit, fmt_len=12, offset=25)
            call plp_dist%separator()
            call plp_dist%print_header()
            do mi = 1, group_size
               call plp_dist%begin_row()
               write (row_label, '(i0)') mi
               call plp_dist%add(trim(row_label))
               do mj = 1, group_size
                  d_ij = norm2(self%xyz(:, igroup_start + mi - 1) &
                               - self%xyz(:, igroup_start + mj - 1))
                  call plp_dist%add(d_ij, fmt='es12.4')
               end do
               call plp_dist%end_row()
            end do
            call plp_dist%separator()
            deallocate (dist_widths, dist_headers)
         end if

         deallocate (phi_local, weights_new, keep)
         igroup_start = igroup_end + 1
      end do

      ! Post-pass summary
      if (post_n_multi_groups > 0) then
         post_avg_kept = real(post_n_multi_kept, wp)/real(post_n_multi_groups, wp)
      else
         post_avg_kept = 0.0_wp
         post_nbranch_max = 0
      end if

      if (self%verbosity >= 2) then
         call pp%push('Branch weighting (final):')
         call pp%kv('Branched anchors', post_n_multi_groups, '(still multi-branch)')
         call pp%kv('Collapsed to single', post_n_collapsed_groups, '(no longer branched)')
         if (post_n_multi_groups > 0) then
            call pp%kv('Branched points (kept)', post_n_multi_kept)
         end if
         call pp%kv('Branched points (dropped)', post_n_dropped_points)
         if (post_n_multi_groups > 0) then
            call pp%kv('Branches/anchor (avg)', post_avg_kept)
            call pp%kv('Branches/anchor (max)', post_nbranch_max)
         end if
         call pp%pop()
      end if

   end subroutine compute_branch_weights

end submodule moist_cavity_drop_filter
