module moist_cavity_drop_types
   use mctc_env_accuracy, only: wp
   use moist_utils_mem, only: grow_array

   implicit none

   public :: projection_workspace_type
   public :: projection_buffer_type

   !> Reusable work arrays for per-anchor projection results
   !>
   !> A workspace is intended to be thread-local and reused across many
   !> projector calls to avoid repeated allocation/deallocation in hot loops
   type :: projection_workspace_type
      !> Number of valid branch entries currently stored.
      integer :: n_points = 0
      !> Allocated branch capacity.
      integer :: capacity = 0

      !> Projected branch coordinates (3, capacity).
      real(wp), allocatable :: points(:, :)
      !> Projected branch normal vectors (3, capacity).
      real(wp), allocatable :: normals(:, :)
      !> Anchor-to-projected displacement norms.
      real(wp), allocatable :: rho(:)
      !> Projection multipliers.
      real(wp), allocatable :: lambda(:)
      !> Objective values (used for branch weighting).
      real(wp), allocatable :: phi(:)
      !> Per-branch weights.
      real(wp), allocatable :: branch_weights(:)
      !> Per-branch convergence flags.
      logical, allocatable :: converged(:)
   contains
      procedure :: init => projection_workspace_init
      procedure :: clear => projection_workspace_clear
      procedure :: reserve => projection_workspace_reserve
      procedure :: set_single => projection_workspace_set_single
      procedure :: size => projection_workspace_size
      procedure :: destroy => projection_workspace_destroy
      final :: finalize_projection_workspace
   end type projection_workspace_type

   !> Thread-local append-only storage for projected DROP grid points
   !>
   !> Stores one entry per projected branch. The buffer keeps capacity separate from logical size (`n_used`)
   !> so parallel projection loops can append efficiently and merge later
   type :: projection_buffer_type
      !> Number of valid entries currently stored in the buffer.
      integer :: n_used = 0
      !> Allocated storage capacity.
      integer :: capacity = 0

      !> Projected grid point positions (3, capacity).
      real(wp), allocatable :: xyz(:, :)
      !> Anchor positions used as projection sources (3, capacity).
      real(wp), allocatable :: anchorxyz(:, :)
      !> Surface normal vectors at projected points (3, capacity).
      real(wp), allocatable :: normal0(:, :)

      !> Projected Lebedev weights.
      real(wp), allocatable :: wleb(:)
      !> Anchor Lebedev weights.
      real(wp), allocatable :: anchor_wleb0(:)
      !> Lagrange multipliers from projection solve.
      real(wp), allocatable :: lambda0(:)
      !> iSwiG anchor switching values.
      real(wp), allocatable :: iswig_f0(:)
      !> Combined switching values.
      real(wp), allocatable :: f(:)
      !> Anchor Gaussian widths.
      real(wp), allocatable :: anchor_xi0(:)
      !> Anchor-to-projected displacement norm.
      real(wp), allocatable :: rho(:)
      !> Branch weights for multi-solution projections.
      real(wp), allocatable :: wbranch(:)
      !> Per-branch objective values at the projected points.
      real(wp), allocatable :: phi0(:)

      !> Owning atom index per point.
      integer, allocatable :: owner(:)
      !> Branch index in anchor group.
      integer, allocatable :: branch(:)
      !> Anchor group id.
      integer, allocatable :: anchor_id(:)
      !> Number of branches for the originating anchor.
      integer, allocatable :: branch_count(:)

      !> Per-point projection convergence flag.
      logical, allocatable :: converged(:)
   contains
      procedure :: init => projection_buffer_init
      procedure :: clear => projection_buffer_clear
      procedure :: reserve => projection_buffer_reserve
      procedure :: append_branches => projection_buffer_append_branches
      procedure :: add_workspace => projection_buffer_add_workspace
      procedure :: size => projection_buffer_size
      procedure :: destroy => projection_buffer_destroy
      final :: finalize_projection_buffer
   end type projection_buffer_type

contains

   !> Compute next storage capacity using +10% growth steps
   !>
   !> @param[in] current_capacity Current allocated capacity
   !> @param[in] required_capacity Minimum required capacity
   !> @param[out] new_capacity Capacity satisfying required capacity
   pure integer function projection_buffer_grow_capacity(current_capacity, required_capacity) result(new_capacity)
      integer, intent(in) :: current_capacity
      integer, intent(in) :: required_capacity

      integer :: increment

      new_capacity = max(1, current_capacity)
      do while (new_capacity < required_capacity)
         increment = max(1, ceiling(0.1_wp*real(new_capacity, wp)))
         new_capacity = new_capacity + increment
      end do
   end function projection_buffer_grow_capacity

   !> Initialize projection workspace and optionally preallocate storage
   !>
   !> @param[inout] self Projection workspace instance
   !> @param[in]    initial_capacity Optional initial branch capacity
   subroutine projection_workspace_init(self, initial_capacity)
      class(projection_workspace_type), intent(inout) :: self
      integer, intent(in), optional :: initial_capacity

      call self%destroy()
      if (present(initial_capacity)) then
         if (initial_capacity > 0) call self%reserve(initial_capacity)
      end if
   end subroutine projection_workspace_init

   !> Reset workspace logical size while keeping allocated memory
   !>
   !> @param[inout] self Projection workspace instance
   subroutine projection_workspace_clear(self)
      class(projection_workspace_type), intent(inout) :: self
      self%n_points = 0
   end subroutine projection_workspace_clear

   !> Ensure workspace can store at least `required_capacity` branches
   !>
   !> @param[inout] self Projection workspace instance
   !> @param[in]    required_capacity Minimum required capacity
   subroutine projection_workspace_reserve(self, required_capacity)
      class(projection_workspace_type), intent(inout) :: self
      integer, intent(in) :: required_capacity
      integer :: new_capacity

      if (required_capacity <= self%capacity) return

      new_capacity = projection_buffer_grow_capacity(self%capacity, required_capacity)

      call grow_array(self%points, 3, new_capacity, fill_value=0.0_wp)
      call grow_array(self%normals, 3, new_capacity, fill_value=0.0_wp)
      call grow_array(self%rho, new_capacity, fill_value=0.0_wp)
      call grow_array(self%lambda, new_capacity, fill_value=0.0_wp)
      call grow_array(self%phi, new_capacity, fill_value=0.0_wp)
      call grow_array(self%branch_weights, new_capacity, fill_value=1.0_wp)
      call grow_array(self%converged, new_capacity, fill_value=.false.)

      self%capacity = new_capacity
   end subroutine projection_workspace_reserve

   !> Store a single-branch result in the workspace
   !>
   !> @param[inout] self Projection workspace instance
   !> @param[in]    point Branch point coordinates
   !> @param[in]    converged_flag Convergence flag for the branch
   subroutine projection_workspace_set_single(self, point, converged_flag)
      class(projection_workspace_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      logical, intent(in) :: converged_flag

      call self%reserve(1)
      self%n_points = 1
      self%points(:, 1) = point
      self%normals(:, 1) = 0.0_wp
      self%rho(1) = 0.0_wp
      self%lambda(1) = 0.0_wp
      self%phi(1) = 0.0_wp
      self%branch_weights(1) = 1.0_wp
      self%converged(1) = converged_flag
   end subroutine projection_workspace_set_single

   !> Return number of valid branch entries currently stored.
   !>
   !> @param[in] self Projection workspace instance
   !> @param[out] n Number of valid branch entries
   pure integer function projection_workspace_size(self) result(n)
      class(projection_workspace_type), intent(in) :: self
      n = self%n_points
   end function projection_workspace_size

   !> Release all workspace allocations
   !>
   !> @param[inout] self Projection workspace instance
   subroutine projection_workspace_destroy(self)
      class(projection_workspace_type), intent(inout) :: self

      self%n_points = 0
      self%capacity = 0

      if (allocated(self%points)) deallocate (self%points)
      if (allocated(self%normals)) deallocate (self%normals)
      if (allocated(self%rho)) deallocate (self%rho)
      if (allocated(self%lambda)) deallocate (self%lambda)
      if (allocated(self%phi)) deallocate (self%phi)
      if (allocated(self%branch_weights)) deallocate (self%branch_weights)
      if (allocated(self%converged)) deallocate (self%converged)
   end subroutine projection_workspace_destroy

   !> Finalizer for projection_workspace_type.
   !>
   !> @param[inout] self Projection workspace instance
   subroutine finalize_projection_workspace(self)
      type(projection_workspace_type), intent(inout) :: self
      call self%destroy()
   end subroutine finalize_projection_workspace

   !> Initialize projection buffer and optionally preallocate storage
   !>
   !> @param[inout] self Projection buffer instance
   !> @param[in]    initial_capacity Optional initial capacity
   subroutine projection_buffer_init(self, initial_capacity)
      class(projection_buffer_type), intent(inout) :: self
      integer, intent(in), optional :: initial_capacity

      call self%destroy()
      if (present(initial_capacity)) then
         if (initial_capacity > 0) call self%reserve(initial_capacity)
      end if
   end subroutine projection_buffer_init

   !> Reset logical size while keeping allocated capacity for reuse
   !>
   !> @param[inout] self Projection buffer instance
   subroutine projection_buffer_clear(self)
      class(projection_buffer_type), intent(inout) :: self
      self%n_used = 0
   end subroutine projection_buffer_clear

   !> Ensure projection buffer has storage for at least `required_capacity`
   !>
   !> @param[inout] self Projection buffer instance
   !> @param[in]    required_capacity Minimum capacity to guarantee
   subroutine projection_buffer_reserve(self, required_capacity)
      class(projection_buffer_type), intent(inout) :: self
      integer, intent(in) :: required_capacity

      integer :: new_capacity

      if (required_capacity <= self%capacity) return

      new_capacity = projection_buffer_grow_capacity(self%capacity, required_capacity)

      call grow_array(self%xyz, 3, new_capacity, fill_value=0.0_wp)
      call grow_array(self%anchorxyz, 3, new_capacity, fill_value=0.0_wp)
      call grow_array(self%normal0, 3, new_capacity, fill_value=0.0_wp)

      call grow_array(self%wleb, new_capacity, fill_value=0.0_wp)
      call grow_array(self%anchor_wleb0, new_capacity, fill_value=0.0_wp)
      call grow_array(self%lambda0, new_capacity, fill_value=0.0_wp)
      call grow_array(self%iswig_f0, new_capacity, fill_value=1.0_wp)
      call grow_array(self%f, new_capacity, fill_value=1.0_wp)
      call grow_array(self%anchor_xi0, new_capacity, fill_value=0.0_wp)
      call grow_array(self%rho, new_capacity, fill_value=0.0_wp)
      call grow_array(self%wbranch, new_capacity, fill_value=1.0_wp)
      call grow_array(self%phi0, new_capacity, fill_value=0.0_wp)

      call grow_array(self%owner, new_capacity, fill_value=0)
      call grow_array(self%branch, new_capacity, fill_value=1)
      call grow_array(self%anchor_id, new_capacity, fill_value=0)
      call grow_array(self%branch_count, new_capacity, fill_value=1)

      call grow_array(self%converged, new_capacity, fill_value=.false.)

      self%capacity = new_capacity
   end subroutine projection_buffer_reserve

   !> Append all projection branches associated with one anchor point
   !>
   !> @param[inout] self Projection buffer instance
   !> @param[in]    anchor_xyz Anchor point coordinates
   !> @param[in]    owner Owning atom index
   !> @param[in]    wleb Anchor Lebedev weight
   !> @param[in]    anchor_wleb0 Raw anchor Lebedev weight
   !> @param[in]    iswig_f0 iSwiG anchor switching value
   !> @param[in]    f Combined switching value
   !> @param[in]    anchor_xi0 Anchor Gaussian width
   !> @param[in]    anchor_id Anchor group id
   !> @param[in]    proj_points Projected branch coordinates
   !> @param[in]    proj_rho Anchor-to-projection distances
   !> @param[in]    proj_lambda Projection multipliers
   !> @param[in]    proj_normals Projected branch normal vectors
   !> @param[in]    proj_converged Branch convergence flags
   !> @param[in]    branch_weights Branch weights
   !> @param[out]   ok Append success flag (optional)
   subroutine projection_buffer_append_branches(self, anchor_xyz, owner, wleb, anchor_wleb0, &
                                                iswig_f0, f, anchor_xi0, anchor_id, proj_points, proj_rho, proj_lambda, &
                                                proj_normals, proj_converged, branch_weights, proj_phi, ok)
      class(projection_buffer_type), intent(inout) :: self
      real(wp), intent(in) :: anchor_xyz(3)
      integer, intent(in) :: owner
      real(wp), intent(in) :: wleb
      real(wp), intent(in) :: anchor_wleb0
      real(wp), intent(in) :: iswig_f0
      real(wp), intent(in) :: f
      real(wp), intent(in) :: anchor_xi0
      integer, intent(in) :: anchor_id
      real(wp), intent(in) :: proj_points(:, :)
      real(wp), intent(in) :: proj_rho(:)
      real(wp), intent(in) :: proj_lambda(:)
      real(wp), intent(in) :: proj_normals(:, :)
      logical, intent(in) :: proj_converged(:)
      real(wp), intent(in) :: branch_weights(:)
      real(wp), intent(in) :: proj_phi(:)
      logical, intent(out), optional :: ok

      integer :: ib, idx, n_branch, first_idx, last_idx
      logical :: valid

      n_branch = size(proj_rho)
      valid = n_branch > 0
      valid = valid .and. size(proj_points, dim=1) == 3
      valid = valid .and. size(proj_points, dim=2) == n_branch
      valid = valid .and. size(proj_lambda) == n_branch
      valid = valid .and. size(proj_normals, dim=1) == 3
      valid = valid .and. size(proj_normals, dim=2) == n_branch
      valid = valid .and. size(proj_converged) == n_branch
      valid = valid .and. size(branch_weights) == n_branch
      valid = valid .and. size(proj_phi) == n_branch
      if (.not. valid) then
         if (present(ok)) ok = .false.
         return
      end if

      first_idx = self%n_used + 1
      last_idx = self%n_used + n_branch
      call self%reserve(last_idx)

      do ib = 1, n_branch
         idx = first_idx + ib - 1
         self%anchorxyz(:, idx) = anchor_xyz
         self%owner(idx) = owner
         self%wleb(idx) = wleb
         self%anchor_wleb0(idx) = anchor_wleb0
         self%iswig_f0(idx) = iswig_f0
         self%f(idx) = f
         self%anchor_xi0(idx) = anchor_xi0
         self%xyz(:, idx) = proj_points(:, ib)
         self%rho(idx) = proj_rho(ib)
         self%lambda0(idx) = proj_lambda(ib)
         self%normal0(:, idx) = proj_normals(:, ib)
         self%converged(idx) = proj_converged(ib)
         self%branch(idx) = ib
         self%anchor_id(idx) = anchor_id
         self%branch_count(idx) = n_branch
         self%wbranch(idx) = branch_weights(ib)
         self%phi0(idx) = proj_phi(ib)
      end do

      self%n_used = last_idx
      if (present(ok)) ok = .true.
   end subroutine projection_buffer_append_branches

   !> Append branches from a projection workspace
   !>
   !> @param[inout] self Projection buffer instance
   !> @param[in]    work Projection workspace containing branch data
   !> @param[in]    anchor_xyz Anchor point coordinates
   !> @param[in]    owner Owning atom index
   !> @param[in]    wleb Anchor Lebedev weight
   !> @param[in]    anchor_wleb0 Raw anchor Lebedev weight
   !> @param[in]    iswig_f0 iSwiG anchor switching value
   !> @param[in]    f Combined switching value
   !> @param[in]    anchor_xi0 Anchor Gaussian width
   !> @param[in]    anchor_id Anchor group id
   !> @param[out]   ok Append success flag (optional)
   subroutine projection_buffer_add_workspace(self, work, anchor_xyz, owner, wleb, anchor_wleb0, &
                                              iswig_f0, f, anchor_xi0, anchor_id, ok)
      class(projection_buffer_type), intent(inout) :: self
      type(projection_workspace_type), intent(in) :: work
      real(wp), intent(in) :: anchor_xyz(3)
      integer, intent(in) :: owner
      real(wp), intent(in) :: wleb
      real(wp), intent(in) :: anchor_wleb0
      real(wp), intent(in) :: iswig_f0
      real(wp), intent(in) :: f
      real(wp), intent(in) :: anchor_xi0
      integer, intent(in) :: anchor_id
      logical, intent(out), optional :: ok
      integer :: n_branch

      n_branch = work%size()
      if (n_branch <= 0) then
         if (present(ok)) ok = .false.
         return
      end if

      call self%append_branches( &
         anchor_xyz=anchor_xyz, &
         owner=owner, &
         wleb=wleb, &
         anchor_wleb0=anchor_wleb0, &
         iswig_f0=iswig_f0, &
         f=f, &
         anchor_xi0=anchor_xi0, &
         anchor_id=anchor_id, &
         proj_points=work%points(:, 1:n_branch), &
         proj_rho=work%rho(1:n_branch), &
         proj_lambda=work%lambda(1:n_branch), &
         proj_normals=work%normals(:, 1:n_branch), &
         proj_converged=work%converged(1:n_branch), &
         branch_weights=work%branch_weights(1:n_branch), &
         proj_phi=work%phi(1:n_branch), &
         ok=ok)
   end subroutine projection_buffer_add_workspace

   !> Return current number of stored projection entries
   !>
   !> @param[in] self Projection buffer instance
   !> @param[out] n Number of used entries
   pure integer function projection_buffer_size(self) result(n)
      class(projection_buffer_type), intent(in) :: self
      n = self%n_used
   end function projection_buffer_size

   !> Destroy buffer content and release all allocated memory
   !>
   !> @param[inout] self Projection buffer instance
   subroutine projection_buffer_destroy(self)
      class(projection_buffer_type), intent(inout) :: self

      self%n_used = 0
      self%capacity = 0

      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%anchorxyz)) deallocate (self%anchorxyz)
      if (allocated(self%normal0)) deallocate (self%normal0)

      if (allocated(self%wleb)) deallocate (self%wleb)
      if (allocated(self%anchor_wleb0)) deallocate (self%anchor_wleb0)
      if (allocated(self%lambda0)) deallocate (self%lambda0)
      if (allocated(self%iswig_f0)) deallocate (self%iswig_f0)
      if (allocated(self%f)) deallocate (self%f)
      if (allocated(self%anchor_xi0)) deallocate (self%anchor_xi0)
      if (allocated(self%rho)) deallocate (self%rho)
      if (allocated(self%wbranch)) deallocate (self%wbranch)
      if (allocated(self%phi0)) deallocate (self%phi0)

      if (allocated(self%owner)) deallocate (self%owner)
      if (allocated(self%branch)) deallocate (self%branch)
      if (allocated(self%anchor_id)) deallocate (self%anchor_id)
      if (allocated(self%branch_count)) deallocate (self%branch_count)

      if (allocated(self%converged)) deallocate (self%converged)
   end subroutine projection_buffer_destroy

   !> Finalizer for projection_buffer_type.
   !>
   !> @param[inout] self Projection buffer instance
   subroutine finalize_projection_buffer(self)
      type(projection_buffer_type), intent(inout) :: self
      call self%destroy()
   end subroutine finalize_projection_buffer

end module moist_cavity_drop_types
