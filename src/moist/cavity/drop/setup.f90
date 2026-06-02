!> DROP setup and preprocessing routines.
submodule(moist_cavity_drop) moist_cavity_drop_setup
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none

contains

   !* ================================================================================= *!
   !*                             Adjacency list/grid setup                             *!
   !* ================================================================================= *!

   !> Set up neighbour list for spheres and build the atom-screening cell grid
   !>
   !> @param[inout] self Cavity instance with molecular structure and parameters
   module subroutine setup_mol_cell_grid(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      real(wp), allocatable :: r_eff(:)
      type(prettyprinter) :: pp

      allocate (r_eff(self%nsph))
      do i = 1, self%nsph
         r_eff(i) = self%radii(i) + self%lsf_model%neighbor_cutoff(self%radii(i))
      end do

      ! Setup real-space cell grid for screening in LSF evaluation
      call self%mol_cell_grid%build(self%mol%xyz, r_eff, &
                                    full_scan_below=self%param%cell_grid_full_scan_below, &
                                    cell_fraction=self%param%cell_grid_fraction)

      ! Build iSwig atom-atom neighbour list (sorted by distance for early exit)
      call self%iswig%update(self%mol, self%radii, wleb_max=maxval(self%anchor_wleb0))

      if (self%verbosity >= 2) then
         pp = new_prettyprinter(unit=output_unit)
         call pp%blank()
         call pp%push('Molecular cell grid:')
         if (self%mol_cell_grid%full_scan) then
            call pp%kv('Mode', 'full-scan (single cell)')
            call pp%kv('Atoms', self%mol_cell_grid%natoms)
         else
            call pp%kv('Mode', 'spatial binning')
            call pp%kv('Cell fraction', self%mol_cell_grid%cell_fraction)
            call pp%kv('Cell side', self%mol_cell_grid%cell_side, 'bohr')
            call pp%kv('Grid nx', self%mol_cell_grid%nx)
            call pp%kv('Grid ny', self%mol_cell_grid%ny)
            call pp%kv('Grid nz', self%mol_cell_grid%nz)
            call pp%kv('Total cells', self%mol_cell_grid%ncells)
            if (self%mol_cell_grid%ncells > 0) then
               call pp%kv('Maximum atoms per cell', maxval(self%mol_cell_grid%cell_nnl))
               call pp%kv('Minimum atoms per cell', minval(self%mol_cell_grid%cell_nnl))
               call pp%kv('Average atoms per cell', &
                          real(sum(self%mol_cell_grid%cell_nnl), wp) &
                          /real(self%mol_cell_grid%ncells, wp), fmt='(f14.4)')
            end if
            call pp%kv('Total cell entries', size(self%mol_cell_grid%cell_nlat))
            call pp%kv('Entries per atom', &
                       real(size(self%mol_cell_grid%cell_nlat), wp)/real(self%nsph, wp), &
                       fmt='(f14.4)')
         end if
         call pp%pop()
      end if

      deallocate (r_eff)

   end subroutine setup_mol_cell_grid

   !> Build grid-point neighbour list for density computation
   !>
   !> @param[inout] self  Cavity instance (anchorxyz must be filled)
   !> @param[out]   error Error object
   module subroutine setup_grid_adj_list(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: cutoff
      type(prettyprinter) :: pp

      cutoff = self%param%adj_list_grid_cutoff

      call self%grid_adj_list%init(cutoff=cutoff)
      call self%grid_adj_list%update(self%xyz(:, 1:self%ngrid))

      if (self%verbosity >= 2) then
         pp = new_prettyprinter(unit=output_unit)
         call pp%blank()
         call pp%push('Grid neighbour list:')
         call pp%kv('Cutoff distance', cutoff, 'bohr')
         call pp%kv('Number of points', self%ngrid)
         call pp%kv('Maximum neighbours', maxval(self%grid_adj_list%nnl))
         call pp%kv('Minimum neighbours', minval(self%grid_adj_list%nnl))
         call pp%kv('Average neighbours', &
                    real(sum(self%grid_adj_list%nnl), wp)/real(self%ngrid, wp), fmt='(f14.4)')
         call pp%kv('Average neighbours (%)', &
                    real(sum(self%grid_adj_list%nnl), wp)/real(self%ngrid, wp)/real(self%ngrid, wp)*100.0_wp, '%', &
                    fmt='(f14.4)')
         call pp%kv('Total pairs', size(self%grid_adj_list%nlat))
         call pp%kv('Total pairs (%)', &
                    real(size(self%grid_adj_list%nlat), wp)/real(self%ngrid**2, wp)*100.0_wp, '%', fmt='(f14.4)')
         call pp%pop()
      end if

   end subroutine setup_grid_adj_list

   !* ================================================================================= *!
   !*                              Initial filling of arrays                             *!
   !* ================================================================================= *!

   !> Initialize grid arrays with Lebedev points on atomic spheres
   !>
   !> Places num_leb points on each atomic sphere using cached Lebedev grid.
   !> Sets xyz, wleb (weights), and owner (atom index) arrays.
   !>
   !> @param[inout] self Cavity instance
   module subroutine fill_arrays(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      integer :: i, ii, jj

      ! Keep anchor-grid size consistent even after branched projection resized nmax.
      self%nmax = self%nsph*self%param%num_leb
      self%ngrid = self%nmax

      ! Allocate grid arrays
      if (allocated(self%xyz)) deallocate (self%xyz)
      allocate (self%xyz(3, self%nmax), source=0.0_wp)
      if (allocated(self%lambda0)) deallocate (self%lambda0)
      allocate (self%lambda0(self%nmax), source=0.0_wp)
      if (allocated(self%anchorxyz)) deallocate (self%anchorxyz)
      allocate (self%anchorxyz(3, self%nmax), source=0.0_wp)
      if (allocated(self%wleb)) deallocate (self%wleb)
      allocate (self%wleb(self%nmax), source=0.0_wp)
      if (allocated(self%owner)) deallocate (self%owner)
      allocate (self%owner(self%nmax), source=0)
      if (allocated(self%numbering)) deallocate (self%numbering)
      if (allocated(self%anchor_wleb0)) deallocate (self%anchor_wleb0)
      allocate (self%anchor_wleb0(self%nmax), source=0.0_wp)
      if (allocated(self%branch)) deallocate (self%branch)
      allocate (self%branch(self%nmax), source=1)
      if (allocated(self%anchor_id)) deallocate (self%anchor_id)
      allocate (self%anchor_id(self%nmax), source=0)
      if (allocated(self%branch_count)) deallocate (self%branch_count)
      allocate (self%branch_count(self%nmax), source=1)
      if (allocated(self%wbranch)) deallocate (self%wbranch)
      allocate (self%wbranch(self%nmax), source=1.0_wp)
      if (allocated(self%phi0)) deallocate (self%phi0)
      if (allocated(self%normal0)) deallocate (self%normal0)
      allocate (self%normal0(3, self%nmax), source=0.0_wp)
      if (allocated(self%rho)) deallocate (self%rho)
      allocate (self%rho(self%nmax), source=0.0_wp)
      if (allocated(self%r_iI0)) deallocate (self%r_iI0)
      allocate (self%r_iI0(self%nmax), source=0.0_wp)
      if (allocated(self%converged)) deallocate (self%converged)
      allocate (self%converged(self%nmax), source=.false.)
      if (allocated(self%cpjac_scal0)) deallocate (self%cpjac_scal0)
      if (allocated(self%w_f0)) deallocate (self%w_f0)

      ii = 0
      do i = 1, self%nsph
         do jj = 1, self%param%num_leb
            ii = ii + 1

            ! Construct raw Lebedev weight from ang_weight(jj)
            self%wleb(ii) = self%ang_weight(jj)*(4.0_wp*pi)

            ! Cartesian location of point on sphere i:
            self%xyz(1, ii) = self%mol%xyz(1, i) + self%radii(i)*self%ang_grid(1, jj)
            self%xyz(2, ii) = self%mol%xyz(2, i) + self%radii(i)*self%ang_grid(2, jj)
            self%xyz(3, ii) = self%mol%xyz(3, i) + self%radii(i)*self%ang_grid(3, jj)

            ! owner-atom index:
            self%owner(ii) = i
            self%anchor_id(ii) = ii
         end do
      end do

      self%anchorxyz = self%xyz
      self%anchor_wleb0 = self%wleb

   end subroutine fill_arrays

   !* ================================================================================= *!
   !*                                     Switching                                     *!
   !* ================================================================================= *!

   module subroutine compute_switching_function(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      type(prettyprinter) :: pp
      real(wp) :: lsf0_gradnorm
      real(wp) :: perc_iswig, perc_anchor

      integer :: i
      integer :: cg_start, cg_n
      integer :: n_iswig_removed, n_anchor_additional
      character(64) :: stat_line

      if (allocated(self%anchor_xi0)) deallocate (self%anchor_xi0)
      allocate (self%anchor_xi0(self%nmax), source=0.0_wp)
      if (allocated(self%iswig_f0)) deallocate (self%iswig_f0)
      allocate (self%iswig_f0(self%nmax), source=1.0_wp)
      if (allocated(self%f)) deallocate (self%f)
      allocate (self%f(self%nmax), source=1.0_wp)
      !$omp parallel default(shared) private(i, lsf, lsf0_gradnorm, cg_start, cg_n)
      allocate (lsf, source=self%lsf_model)

      !$omp do schedule(static)
      do i = 1, self%nmax

         ! iswig (uses built-in sorted neighbor list with early exit)
         self%anchor_xi0(i) = self%iswig%xi0(self%owner(i), self%anchor_wleb0(i))
         self%iswig_f0(i) = self%iswig%swi0(self%anchorxyz(:, i), self%owner(i), self%anchor_xi0(i))

      end do
      !$omp end do
      !$omp end parallel

      !> Total switching function
      self%f = self%iswig_f0

      if (self%verbosity > 1) then
         n_iswig_removed = count(self%iswig_f0 < self%param%wleb_cut)
         n_anchor_additional = count((self%iswig_f0 >= self%param%wleb_cut) .and. &
                                     (self%f < self%param%wleb_cut))

         if (self%nmax > 0) then
            perc_iswig = real(n_iswig_removed, wp)/real(self%nmax, wp)*100.0_wp
            perc_anchor = real(n_anchor_additional, wp)/real(self%nmax, wp)*100.0_wp
         else
            perc_iswig = 0.0_wp
            perc_anchor = 0.0_wp
         end if

         pp = new_prettyprinter(unit=output_unit)
         call pp%blank()
         call pp%push('Pre-filter:')
         call pp%kv('Points before removal', self%nmax)
         call pp%kv('iSwig removed', n_iswig_removed)
         call pp%kv('anchor_f additional', n_anchor_additional)
         call pp%kv('Points to project', self%nmax - n_iswig_removed - n_anchor_additional)
         call pp%pop()
      end if

   end subroutine compute_switching_function

end submodule moist_cavity_drop_setup
