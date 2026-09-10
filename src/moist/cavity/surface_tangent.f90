!> Directional tangents of the common quantities of a discretized cavity surface
!>
!> The forward-mode dual of [[cavity_surface_adjoint_type]]: where the adjoint
!> accumulator carries one weight per surface observable, this carries the
!> directional derivative of every observable along a batch of nuclear
!> directions, `d(Gamma)/dR . v` for each direction `v`. A cavity fills it
!> through `get_surface_tangent`; a model component reads it to differentiate
!> its surface weights along the same directions, which is the second-order
!> surface channel of the nuclear Hessian.
!>
!> Every channel is `(ngrid, ndir)` or `(3, ngrid, ndir)`; column `idir` holds
!> the tangent along direction `idir` of the batch. The curvature channels are
!> meaningful only when `have_curvature` is set: a cavity that was not asked for
!> curvature tangents leaves them allocated and zero, and a reader that needs
!> them must check the flag rather than trust the zeros.
module moist_cavity_surface_tangent
   use mctc_env, only: wp

   implicit none (type, external)
   private

   public :: cavity_surface_tangent_type

   !> Directional tangents of the common quantities of a discretized cavity surface
   type :: cavity_surface_tangent_type
      !> Tangents of the Gaussian widths xi_i (ngrid, ndir)
      real(wp), allocatable :: d_xi(:, :)
      !> Tangents of the switching factors f_i (ngrid, ndir)
      real(wp), allocatable :: d_f(:, :)
      !> Tangents of the surface areas a_i (ngrid, ndir)
      real(wp), allocatable :: d_a(:, :)
      !> Tangents of the integration weights w_i (ngrid, ndir)
      real(wp), allocatable :: d_w(:, :)
      !> Tangents of the surface positions r_i (3, ngrid, ndir)
      real(wp), allocatable :: d_xyz(:, :, :)
      !> Tangents of the outward normals n_i (3, ngrid, ndir)
      real(wp), allocatable :: d_n(:, :, :)
      !> Tangents of the first principal curvature k1_i (ngrid, ndir)
      real(wp), allocatable :: d_k1(:, :)
      !> Tangents of the second principal curvature k2_i (ngrid, ndir)
      real(wp), allocatable :: d_k2(:, :)
      !> Whether the curvature channels were requested and are meaningful
      logical :: have_curvature = .false.
   contains
      !> Allocate every channel to zero for a grid and a direction batch
      procedure :: init => init_surface_tangent
      !> Report whether every channel is allocated with consistent shapes
      procedure :: is_initialized => surface_tangent_is_initialized
      !> Number of directions of the batch, zero when uninitialized
      procedure :: ndir => surface_tangent_ndir
   end type cavity_surface_tangent_type

contains

   !> Allocate every channel to zero for a grid and a direction batch
   !>
   !> @param[inout] self           Surface tangent
   !> @param[in]    ngrid          Number of surface grid points
   !> @param[in]    ndir           Number of nuclear directions of the batch
   !> @param[in]    want_curvature Whether the curvature channels are to be filled
   subroutine init_surface_tangent(self, ngrid, ndir, want_curvature)
      !> Surface tangent
      class(cavity_surface_tangent_type), intent(inout) :: self
      !> Number of surface grid points
      integer, intent(in) :: ngrid
      !> Number of nuclear directions of the batch
      integer, intent(in) :: ndir
      !> Whether the curvature channels are to be filled
      logical, intent(in) :: want_curvature

      if (allocated(self%d_xi)) deallocate (self%d_xi)
      if (allocated(self%d_f)) deallocate (self%d_f)
      if (allocated(self%d_a)) deallocate (self%d_a)
      if (allocated(self%d_w)) deallocate (self%d_w)
      if (allocated(self%d_xyz)) deallocate (self%d_xyz)
      if (allocated(self%d_n)) deallocate (self%d_n)
      if (allocated(self%d_k1)) deallocate (self%d_k1)
      if (allocated(self%d_k2)) deallocate (self%d_k2)

      allocate (self%d_xi(ngrid, ndir), source=0.0_wp)
      allocate (self%d_f(ngrid, ndir), source=0.0_wp)
      allocate (self%d_a(ngrid, ndir), source=0.0_wp)
      allocate (self%d_w(ngrid, ndir), source=0.0_wp)
      allocate (self%d_xyz(3, ngrid, ndir), source=0.0_wp)
      allocate (self%d_n(3, ngrid, ndir), source=0.0_wp)
      allocate (self%d_k1(ngrid, ndir), source=0.0_wp)
      allocate (self%d_k2(ngrid, ndir), source=0.0_wp)
      self%have_curvature = want_curvature

   end subroutine init_surface_tangent

   !> Check whether every channel has been allocated consistently
   !>
   !> @param[in] self  Surface tangent
   pure logical function surface_tangent_is_initialized(self) result(initialized)
      !> Surface tangent
      class(cavity_surface_tangent_type), intent(in) :: self

      !> Grid and direction extents
      integer :: ngrid, ndir

      initialized = allocated(self%d_xi) .and. allocated(self%d_f) .and. &
                    allocated(self%d_a) .and. allocated(self%d_w) .and. &
                    allocated(self%d_xyz) .and. allocated(self%d_n) .and. &
                    allocated(self%d_k1) .and. allocated(self%d_k2)
      if (.not. initialized) return

      ngrid = size(self%d_xi, 1)
      ndir = size(self%d_xi, 2)
      initialized = all(shape(self%d_f) == [ngrid, ndir]) .and. &
                    all(shape(self%d_a) == [ngrid, ndir]) .and. &
                    all(shape(self%d_w) == [ngrid, ndir]) .and. &
                    all(shape(self%d_xyz) == [3, ngrid, ndir]) .and. &
                    all(shape(self%d_n) == [3, ngrid, ndir]) .and. &
                    all(shape(self%d_k1) == [ngrid, ndir]) .and. &
                    all(shape(self%d_k2) == [ngrid, ndir])

   end function surface_tangent_is_initialized

   !> Number of directions of the batch, zero when uninitialized
   !>
   !> @param[in] self  Surface tangent
   pure integer function surface_tangent_ndir(self) result(ndir)
      !> Surface tangent
      class(cavity_surface_tangent_type), intent(in) :: self

      ndir = 0
      if (allocated(self%d_xi)) ndir = size(self%d_xi, 2)

   end function surface_tangent_ndir

end module moist_cavity_surface_tangent
