!> Owned scratch storage for sequential primal and transpose Krylov solves
module moist_model_moz_solver_jfnk_workspace
   use mctc_env, only: wp
   implicit none(type, external)
   private
   public :: moist_model_moz_solver_workspace_type, exchange_solver_workspace

   !> Owned scratch arrays of a Newton-Krylov solver
   !>
   !> - Allocatable arrays keep their layout when ownership moves between solvers
   !> - Solver objects retain convergence settings, counters and validity flags
   !> - Recycled vectors/images are stored here but invalidated by the solver on transfer
   type :: moist_model_moz_solver_workspace_type
      !> Shape of the real-space fields; zero before allocation
      integer :: ws_nflat = 0, ws_npts = 0, ws_npair = 0
      !> Shared scratch retains capacity for both requested restart dimensions
      logical :: shared = .false.
      !> Newton trial field and correction, each (npts, npair)
      real(wp), allocatable :: ws_gamma_trial(:, :), ws_delta(:, :)
      !> Residuals at the current and trial fields, each (npts, npair)
      real(wp), allocatable :: ws_f_curr(:, :), ws_f_trial(:, :)
      !> Perturbed field and its residual, each (npts, npair)
      !>
      !> - Allocated only when a finite-difference product is actually needed
      real(wp), allocatable :: ws_gamma_pert(:, :), ws_f_pert(:, :)
      !> Flat right-hand side, correction and Arnoldi work vectors, each (nflat)
      real(wp), allocatable :: ws_rhs_flat(:), ws_delta_flat(:), ws_r(:), ws_w(:)
      !> Optional deflation work vector (nflat)
      real(wp), allocatable :: ws_z(:)
      !> Arnoldi basis (nflat, capacity + 1) and Hessenberg matrix (capacity + 1, capacity)
      !>
      !> - Capacity follows the requested restart, not the number of steps actually taken
      real(wp), allocatable :: ws_v(:, :), ws_h(:, :)
      !> Givens cosines and sines (capacity), least-squares right-hand side (capacity + 1) and solution (capacity)
      real(wp), allocatable :: ws_cs(:), ws_sn(:), ws_g(:), ws_y(:)
      !> Projection coefficient work vectors, each (capacity)
      real(wp), allocatable :: ws_proj(:), ws_proj2(:)
      !> Recycled basis and image vectors (nflat, recycle capacity)
      !>
      !> - Validity remains in the solver state
      real(wp), allocatable :: ws_u(:, :), ws_c(:, :)
      !> Recycling Hessenberg copy (capacity + 1, capacity) and projection coefficients (recycle capacity, capacity)
      real(wp), allocatable :: ws_hbar(:, :), ws_bproj(:, :)
   contains
      procedure :: ensure => workspace_ensure
      procedure :: ensure_recycle => workspace_ensure_recycle
      procedure :: ensure_fd => workspace_ensure_fd
   end type moist_model_moz_solver_workspace_type

contains

   !> Reuse compatible scratch or reallocate and first-touch every array
   !>
   !> - Shared scratch: reused with capacity >= m + 1; only recycling storage grows
   !> - Standalone scratch: reused only with capacity m + 1 and recycle capacity kmax
   !> - Reallocation drops every array, including FD and deflation scratch, and keeps `shared`
   !> - True result: storage was (re)allocated and any recycled basis is invalid
   !>
   !> @param[inout] self      owned scratch
   !> @param[in]    npts      grid points per field
   !> @param[in]    npair     site pairs per field
   !> @param[in]    m         restart dimension, positive
   !> @param[in]    kmax      recycled dimension; zero allocates no recycled basis
   !> @param[in]    nthreads  OpenMP team size for the first touch
   function workspace_ensure(self, npts, npair, m, kmax, nthreads) result(reallocated)
      !> Owned scratch
      class(moist_model_moz_solver_workspace_type), intent(inout) :: self
      !> Grid points per field
      integer, intent(in) :: npts
      !> Site pairs per field
      integer, intent(in) :: npair
      !> Restart dimension, positive
      integer, intent(in) :: m
      !> Recycled dimension; zero allocates no recycled basis
      integer, intent(in) :: kmax
      !> OpenMP team size for the first touch
      integer, intent(in) :: nthreads
      !> Whether storage was (re)allocated and any recycled basis is invalid
      logical :: reallocated
      !> Flat field length
      integer :: nflat
      !> Shared flag carried over a reallocation
      logical :: shared

      nflat = npts*npair
      if (self%ws_nflat == nflat .and. self%ws_npts == npts .and. self%ws_npair == npair) then
         if (allocated(self%ws_v)) then
            if (self%shared .and. size(self%ws_v, 2) >= m + 1) then
               reallocated = self%ensure_recycle(kmax, nthreads)
               return
            end if
            reallocated = .false.
            if (size(self%ws_v, 2) == m + 1) then
               if (kmax == 0) then
                  if (.not. allocated(self%ws_u)) return
               else if (allocated(self%ws_u)) then
                  if (size(self%ws_u, 2) == kmax) return
               end if
            end if
         end if
      end if

      shared = self%shared
      call clear_workspace(self)
      self%shared = shared
      call allocate_workspace(self, npts, npair, m, kmax, nthreads)
      reallocated = .true.
   end function workspace_ensure

   !> Grow optional recycling storage without replacing a shared Arnoldi basis
   !>
   !> - Projection rows follow the allocated Arnoldi capacity, not the requested restart
   !> - True result: recycling storage was replaced and any recycled basis is invalid
   !>
   !> @param[inout] self      shape-compatible scratch with an allocated Arnoldi basis
   !> @param[in]    kmax      recycled dimension; zero does not consume storage
   !> @param[in]    nthreads  OpenMP team size for the first touch
   function workspace_ensure_recycle(self, kmax, nthreads) result(reallocated)
      !> Shape-compatible scratch with an allocated Arnoldi basis
      class(moist_model_moz_solver_workspace_type), intent(inout) :: self
      !> Recycled dimension; zero does not consume storage
      integer, intent(in) :: kmax
      !> OpenMP team size for the first touch
      integer, intent(in) :: nthreads
      !> Whether recycling storage was replaced and any recycled basis is invalid
      logical :: reallocated
      !> Loop indices, flat field length and Arnoldi capacity
      integer :: i, j, nflat, m

      reallocated = .false.
      if (kmax == 0) return
      if (allocated(self%ws_u)) then
         if (size(self%ws_u, 2) >= kmax) return
      end if
      nflat = self%ws_nflat
      m = size(self%ws_v, 2) - 1
      if (allocated(self%ws_u)) deallocate (self%ws_u, self%ws_c)
      deallocate (self%ws_bproj)
      allocate (self%ws_u(nflat, kmax), self%ws_c(nflat, kmax))
      allocate (self%ws_bproj(kmax, m), source=0.0_wp)
      !$omp parallel do collapse(2) schedule(static) num_threads(nthreads)
      do j = 1, kmax
         do i = 1, nflat
            self%ws_u(i, j) = 0.0_wp
            self%ws_c(i, j) = 0.0_wp
         end do
      end do
      !$omp end parallel do
      reallocated = .true.
   end function workspace_ensure_recycle

   !> Lazily allocate the finite-difference fields, (npts, npair) each
   !>
   !> - Products overwrite every value, so the fields are not initialized
   !> - Existing fields are kept as they are
   !>
   !> @param[inout] self   owned scratch
   !> @param[in]    npts   grid points per field
   !> @param[in]    npair  site pairs per field
   subroutine workspace_ensure_fd(self, npts, npair)
      !> Owned scratch
      class(moist_model_moz_solver_workspace_type), intent(inout) :: self
      !> Grid points per field
      integer, intent(in) :: npts
      !> Site pairs per field
      integer, intent(in) :: npair

      if (allocated(self%ws_gamma_pert)) return
      allocate (self%ws_gamma_pert(npts, npair), self%ws_f_pert(npts, npair))
   end subroutine workspace_ensure_fd

   !> Deallocate every array and restore default components
   !>
   !> @param[out] ws  scratch reset to its default state
   subroutine clear_workspace(ws)
      !> Scratch reset to its default state
      type(moist_model_moz_solver_workspace_type), intent(out) :: ws
   end subroutine clear_workspace

   !> Allocate all solve arrays of cleared scratch and first-touch them
   !>
   !> - Krylov columns are touched with the same static schedule as later vector sweeps
   !> - FD and deflation scratch stay unallocated until needed
   !>
   !> @param[inout] self      cleared scratch
   !> @param[in]    npts      grid points per field
   !> @param[in]    npair     site pairs per field
   !> @param[in]    m         restart dimension, positive
   !> @param[in]    kmax      recycled dimension; zero allocates no recycled basis
   !> @param[in]    nthreads  OpenMP team size for the first touch
   subroutine allocate_workspace(self, npts, npair, m, kmax, nthreads)
      !> Cleared scratch
      type(moist_model_moz_solver_workspace_type), intent(inout) :: self
      !> Grid points per field
      integer, intent(in) :: npts
      !> Site pairs per field
      integer, intent(in) :: npair
      !> Restart dimension, positive
      integer, intent(in) :: m
      !> Recycled dimension; zero allocates no recycled basis
      integer, intent(in) :: kmax
      !> OpenMP team size for the first touch
      integer, intent(in) :: nthreads
      !> Loop indices and flat field length
      integer :: i, j, nflat

      nflat = npts*npair
      self%ws_nflat = nflat
      self%ws_npts = npts
      self%ws_npair = npair

      allocate (self%ws_gamma_trial(npts, npair))
      allocate (self%ws_delta(npts, npair))
      allocate (self%ws_f_curr(npts, npair))
      allocate (self%ws_f_trial(npts, npair))

      allocate (self%ws_rhs_flat(nflat))
      allocate (self%ws_delta_flat(nflat))

      allocate (self%ws_r(nflat))
      allocate (self%ws_w(nflat))

      allocate (self%ws_v(nflat, m + 1))
      allocate (self%ws_h(m + 1, m))
      allocate (self%ws_cs(m))
      allocate (self%ws_sn(m))
      allocate (self%ws_g(m + 1))
      allocate (self%ws_y(m))
      allocate (self%ws_proj(m))
      allocate (self%ws_proj2(m))
      allocate (self%ws_hbar(m + 1, m), source=0.0_wp)
      allocate (self%ws_bproj(max(1, kmax), m), source=0.0_wp)
      if (kmax > 0) allocate (self%ws_u(nflat, kmax), self%ws_c(nflat, kmax))

      !$omp parallel num_threads(nthreads) private(i, j)
      !$omp do collapse(2) schedule(static)
      do j = 1, npair
         do i = 1, npts
            self%ws_gamma_trial(i, j) = 0.0_wp
            self%ws_delta(i, j) = 0.0_wp
            self%ws_f_curr(i, j) = 0.0_wp
            self%ws_f_trial(i, j) = 0.0_wp
         end do
      end do
      !$omp end do
      !$omp do schedule(static)
      do i = 1, nflat
         self%ws_rhs_flat(i) = 0.0_wp
         self%ws_delta_flat(i) = 0.0_wp
         self%ws_r(i) = 0.0_wp
         self%ws_w(i) = 0.0_wp
      end do
      !$omp end do
      ! Match first touch of every Krylov column to its later vector sweeps
      do j = 1, m + 1
         !$omp do schedule(static)
         do i = 1, nflat
            self%ws_v(i, j) = 0.0_wp
         end do
         !$omp end do
      end do
      do j = 1, kmax
         !$omp do schedule(static)
         do i = 1, nflat
            self%ws_u(i, j) = 0.0_wp
            self%ws_c(i, j) = 0.0_wp
         end do
         !$omp end do
      end do
      !$omp end parallel
      self%ws_h = 0.0_wp
      self%ws_g = 0.0_wp
      self%ws_y = 0.0_wp
   end subroutine allocate_workspace

   !> Swap complete allocations while preserving array layout and shared capacity
   !>
   !> @param[inout] first   owned scratch of one completed solver
   !> @param[inout] second  owned scratch of a distinct solver or solve context
   subroutine exchange_solver_workspace(first, second)
      !> Owned scratch of one completed solver
      type(moist_model_moz_solver_workspace_type), allocatable, intent(inout) :: first
      !> Owned scratch of a distinct solver or solve context
      type(moist_model_moz_solver_workspace_type), allocatable, intent(inout) :: second
      !> Temporary holding the first allocation during the swap
      type(moist_model_moz_solver_workspace_type), allocatable :: held

      call move_alloc(first, held)
      call move_alloc(second, first)
      call move_alloc(held, second)
      if (allocated(first)) first%shared = .true.
      if (allocated(second)) second%shared = .true.
   end subroutine exchange_solver_workspace
end module moist_model_moz_solver_jfnk_workspace
