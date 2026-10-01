!> Volumetric position, quadrature-weight and Gaussian-width adjoints
module moist_math_grid_3d_adjoint
   use mctc_env, only: wp, error_type, fatal_error

   implicit none(type, external)
   private

   public :: volume_adjoint_type

   !> Adjoint weights for volumetric point positions, quadrature weights and widths
   type :: volume_adjoint_type
      !> Weights for point positions r_i (3, ngrid)
      real(wp), allocatable :: w_xyz(:, :)
      !> Weights for integration weights w_i (ngrid)
      real(wp), allocatable :: w_w(:)
      !> Weights for Gaussian widths xi_i (ngrid)
      real(wp), allocatable :: w_xi(:)
   contains
      !> Allocate every channel and initialize it to zero
      procedure :: init => init_domain_adjoint
      !> Reset every allocated channel to zero
      procedure :: zero => zero_domain_adjoint
      !> Add any supplied channels
      procedure :: add_weights => add_domain_weights
      !> Report whether every channel is allocated with consistent shapes
      procedure :: is_initialized => domain_adjoint_is_initialized
      !> Number of grid points the accumulator was initialized for, 0 if none
      procedure :: size => domain_adjoint_size
   end type volume_adjoint_type

contains

   !> Allocate every channel and initialize it to zero
   !>
   !> @param[in,out] self  Accumulator
   !> @param[in]    ngrid Number of grid points
   subroutine init_domain_adjoint(self, ngrid)
      !> Accumulator
      class(volume_adjoint_type), intent(inout) :: self
      !> Number of grid points
      integer, intent(in) :: ngrid

      if (allocated(self%w_xyz)) deallocate (self%w_xyz)
      if (allocated(self%w_w)) deallocate (self%w_w)
      if (allocated(self%w_xi)) deallocate (self%w_xi)

      allocate (self%w_xyz(3, ngrid), source=0.0_wp)
      allocate (self%w_w(ngrid), source=0.0_wp)
      allocate (self%w_xi(ngrid), source=0.0_wp)

   end subroutine init_domain_adjoint

   !> Reset every allocated channel to zero
   !>
   !> @param[in,out] self  Accumulator
   subroutine zero_domain_adjoint(self)
      !> Accumulator
      class(volume_adjoint_type), intent(inout) :: self

      if (allocated(self%w_xyz)) self%w_xyz = 0.0_wp
      if (allocated(self%w_w)) self%w_w = 0.0_wp
      if (allocated(self%w_xi)) self%w_xi = 0.0_wp

   end subroutine zero_domain_adjoint

   !> Add any supplied weights to an initialized accumulator
   !>
   !> @param[in,out] self  Accumulator
   !> @param[out]   error Error handling
   !> @param[in]    w_xyz Optional position weights (3, ngrid)
   !> @param[in]    w_w   Optional integration-weight weights (ngrid)
   !> @param[in]    w_xi  Optional Gaussian-width weights (ngrid)
   subroutine add_domain_weights(self, error, w_xyz, w_w, w_xi)
      !> Accumulator
      class(volume_adjoint_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional position weights
      real(wp), intent(in), optional :: w_xyz(:, :)
      !> Optional scalar weights
      real(wp), intent(in), optional :: w_w(:)
      !> Optional Gaussian-width weights
      real(wp), intent(in), optional :: w_xi(:)

      !> Number of grid points
      integer :: ngrid

      if (.not. domain_adjoint_is_initialized(self)) then
         call fatal_error(error, "add_weights: accumulator is not initialized")
         return
      end if
      ngrid = size(self%w_xi)

      if (present(w_xyz)) then
         if (size(w_xyz, 1) /= 3 .or. size(w_xyz, 2) /= ngrid) then
            call fatal_error(error, "add_weights: xyz weight shape mismatch")
            return
         end if
      end if
      if (present(w_w)) then
         if (size(w_w) /= ngrid) then
            call fatal_error(error, "add_weights: integration weight size mismatch")
            return
         end if
      end if
      if (present(w_xi)) then
         if (size(w_xi) /= ngrid) then
            call fatal_error(error, "add_weights: xi weight size mismatch")
            return
         end if
      end if

      if (present(w_xyz)) self%w_xyz = self%w_xyz + w_xyz
      if (present(w_w)) self%w_w = self%w_w + w_w
      if (present(w_xi)) self%w_xi = self%w_xi + w_xi

   end subroutine add_domain_weights

   !> Check whether every channel has been initialized consistently
   !>
   !> @param[in] self  Accumulator
   pure function domain_adjoint_is_initialized(self) result(initialized)
      !> Accumulator
      class(volume_adjoint_type), intent(in) :: self

      !> Initialization status
      logical :: initialized

      !> Number of grid points
      integer :: ngrid

      initialized = allocated(self%w_xyz) .and. allocated(self%w_w) .and. allocated(self%w_xi)
      if (.not. initialized) return

      ngrid = size(self%w_xi)
      initialized = size(self%w_w) == ngrid .and. &
                    size(self%w_xyz, 1) == 3 .and. size(self%w_xyz, 2) == ngrid

   end function domain_adjoint_is_initialized

   !> Number of grid points the accumulator holds channels for, 0 before `init`
   !>
   !> @param[in] self  Accumulator
   pure function domain_adjoint_size(self) result(ngrid)
      !> Accumulator
      class(volume_adjoint_type), intent(in) :: self
      !> Number of domain points
      integer :: ngrid

      ngrid = 0
      if (allocated(self%w_xi)) ngrid = size(self%w_xi)

   end function domain_adjoint_size

end module moist_math_grid_3d_adjoint
