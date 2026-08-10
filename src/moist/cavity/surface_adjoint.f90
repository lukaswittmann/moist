!> Generic accumulation of adjoints with respect to cavity-surface quantities
module moist_cavity_surface_adjoint
   use mctc_env, only: wp, error_type, fatal_error

   implicit none (type, external)
   private

   public :: cavity_surface_adjoint_type

   !> Adjoint weights for the common quantities of a discretized cavity surface
   type :: cavity_surface_adjoint_type
      !> Weights for Gaussian widths xi_i (ngrid)
      real(wp), allocatable :: w_xi(:)
      !> Weights for switching factors f_i (ngrid)
      real(wp), allocatable :: w_f(:)
      !> Weights for surface areas a_i (ngrid)
      real(wp), allocatable :: w_a(:)
      !> Weights for integration weights w_i (ngrid)
      real(wp), allocatable :: w_w(:)
      !> Weights for surface positions r_i (3, ngrid)
      real(wp), allocatable :: w_xyz(:, :)
      !> Weights for outward normals n_i (3, ngrid)
      real(wp), allocatable :: w_n(:, :)
      !> Weights for the first principal curvature k1_i (ngrid)
      real(wp), allocatable :: w_k1(:)
      !> Weights for the second principal curvature k2_i (ngrid)
      real(wp), allocatable :: w_k2(:)
   contains
      !> Allocate every surface-adjoint channel and initialize it to zero
      procedure :: init => init_surface_adjoint
      !> Reset every allocated surface-adjoint channel to zero
      procedure :: zero => zero_surface_adjoint
      !> Add any supplied surface-adjoint channels
      procedure :: add_surface_weights
      !> Report whether every channel is allocated with consistent shapes
      procedure :: is_initialized => surface_adjoint_is_initialized
   end type cavity_surface_adjoint_type

contains

   !> Allocate every surface-adjoint channel and initialize it to zero
   !>
   !> @param[inout] self  Surface-adjoint accumulator
   !> @param[in]    ngrid Number of surface grid points
   subroutine init_surface_adjoint(self, ngrid)
      !> Surface-adjoint accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: self
      !> Number of surface grid points
      integer, intent(in) :: ngrid

      if (allocated(self%w_xi)) deallocate (self%w_xi)
      if (allocated(self%w_f)) deallocate (self%w_f)
      if (allocated(self%w_a)) deallocate (self%w_a)
      if (allocated(self%w_w)) deallocate (self%w_w)
      if (allocated(self%w_xyz)) deallocate (self%w_xyz)
      if (allocated(self%w_n)) deallocate (self%w_n)
      if (allocated(self%w_k1)) deallocate (self%w_k1)
      if (allocated(self%w_k2)) deallocate (self%w_k2)

      allocate (self%w_xi(ngrid), source=0.0_wp)
      allocate (self%w_f(ngrid), source=0.0_wp)
      allocate (self%w_a(ngrid), source=0.0_wp)
      allocate (self%w_w(ngrid), source=0.0_wp)
      allocate (self%w_xyz(3, ngrid), source=0.0_wp)
      allocate (self%w_n(3, ngrid), source=0.0_wp)
      allocate (self%w_k1(ngrid), source=0.0_wp)
      allocate (self%w_k2(ngrid), source=0.0_wp)

   end subroutine init_surface_adjoint

   !> Reset every allocated surface-adjoint channel to zero
   !>
   !> @param[inout] self  Surface-adjoint accumulator
   subroutine zero_surface_adjoint(self)
      !> Surface-adjoint accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: self

      if (allocated(self%w_xi)) self%w_xi = 0.0_wp
      if (allocated(self%w_f)) self%w_f = 0.0_wp
      if (allocated(self%w_a)) self%w_a = 0.0_wp
      if (allocated(self%w_w)) self%w_w = 0.0_wp
      if (allocated(self%w_xyz)) self%w_xyz = 0.0_wp
      if (allocated(self%w_n)) self%w_n = 0.0_wp
      if (allocated(self%w_k1)) self%w_k1 = 0.0_wp
      if (allocated(self%w_k2)) self%w_k2 = 0.0_wp

   end subroutine zero_surface_adjoint

   !> Add any supplied weights to an initialized surface-adjoint accumulator
   !>
   !> @param[inout] self  Surface-adjoint accumulator
   !> @param[out]   error Error handling
   !> @param[in]    w_xi  Optional Gaussian-width weights (ngrid)
   !> @param[in]    w_f   Optional switching-factor weights (ngrid)
   !> @param[in]    w_xyz Optional surface-position weights (3, ngrid)
   !> @param[in]    w_n   Optional outward-normal weights (3, ngrid)
   !> @param[in]    w_k1  Optional first-principal-curvature weights (ngrid)
   !> @param[in]    w_k2  Optional second-principal-curvature weights (ngrid)
   !> @param[in]    w_a   Optional surface-area weights (ngrid)
   !> @param[in]    w_w   Optional integration-weight weights (ngrid)
   subroutine add_surface_weights(self, error, w_xi, w_f, w_xyz, w_n, w_k1, w_k2, w_a, w_w)
      !> Surface-adjoint accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional scalar surface weights
      real(wp), intent(in), optional :: w_xi(:), w_f(:), w_k1(:), w_k2(:), w_a(:), w_w(:)
      !> Optional vector surface weights
      real(wp), intent(in), optional :: w_xyz(:, :), w_n(:, :)

      !> Number of surface grid points
      integer :: ngrid

      if (.not. surface_adjoint_is_initialized(self)) then
         call fatal_error(error, "add_surface_weights: accumulator is not initialized")
         return
      end if
      ngrid = size(self%w_xi)

      if (present(w_xi)) then
         if (size(w_xi) /= ngrid) then
            call fatal_error(error, "add_surface_weights: xi weight size mismatch")
            return
         end if
      end if
      if (present(w_f)) then
         if (size(w_f) /= ngrid) then
            call fatal_error(error, "add_surface_weights: f weight size mismatch")
            return
         end if
      end if
      if (present(w_xyz)) then
         if (size(w_xyz, 1) /= 3 .or. size(w_xyz, 2) /= ngrid) then
            call fatal_error(error, "add_surface_weights: xyz weight shape mismatch")
            return
         end if
      end if
      if (present(w_n)) then
         if (size(w_n, 1) /= 3 .or. size(w_n, 2) /= ngrid) then
            call fatal_error(error, "add_surface_weights: normal weight shape mismatch")
            return
         end if
      end if
      if (present(w_k1)) then
         if (size(w_k1) /= ngrid) then
            call fatal_error(error, "add_surface_weights: k1 weight size mismatch")
            return
         end if
      end if
      if (present(w_k2)) then
         if (size(w_k2) /= ngrid) then
            call fatal_error(error, "add_surface_weights: k2 weight size mismatch")
            return
         end if
      end if
      if (present(w_a)) then
         if (size(w_a) /= ngrid) then
            call fatal_error(error, "add_surface_weights: area weight size mismatch")
            return
         end if
      end if
      if (present(w_w)) then
         if (size(w_w) /= ngrid) then
            call fatal_error(error, "add_surface_weights: integration weight size mismatch")
            return
         end if
      end if

      if (present(w_xi)) self%w_xi = self%w_xi + w_xi
      if (present(w_f)) self%w_f = self%w_f + w_f
      if (present(w_xyz)) self%w_xyz = self%w_xyz + w_xyz
      if (present(w_n)) self%w_n = self%w_n + w_n
      if (present(w_k1)) self%w_k1 = self%w_k1 + w_k1
      if (present(w_k2)) self%w_k2 = self%w_k2 + w_k2
      if (present(w_a)) self%w_a = self%w_a + w_a
      if (present(w_w)) self%w_w = self%w_w + w_w

   end subroutine add_surface_weights

   !> Check whether every surface-adjoint channel has been initialized consistently
   !>
   !> @param[in] self  Surface-adjoint accumulator
   pure logical function surface_adjoint_is_initialized(self) result(initialized)
      !> Surface-adjoint accumulator
      class(cavity_surface_adjoint_type), intent(in) :: self

      !> Number of surface grid points
      integer :: ngrid

      initialized = allocated(self%w_xi) .and. allocated(self%w_f) .and. &
                    allocated(self%w_a) .and. allocated(self%w_w) .and. &
                    allocated(self%w_xyz) .and. allocated(self%w_n) .and. &
                    allocated(self%w_k1) .and. allocated(self%w_k2)
      if (.not. initialized) return

      ngrid = size(self%w_xi)
      initialized = size(self%w_f) == ngrid .and. &
                    size(self%w_a) == ngrid .and. size(self%w_w) == ngrid .and. &
                    size(self%w_xyz, 1) == 3 .and. size(self%w_xyz, 2) == ngrid .and. &
                    size(self%w_n, 1) == 3 .and. size(self%w_n, 2) == ngrid .and. &
                    size(self%w_k1) == ngrid .and. size(self%w_k2) == ngrid

   end function surface_adjoint_is_initialized

end module moist_cavity_surface_adjoint
