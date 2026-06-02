!> Closest-point projection objective phi for DROP cavity grid projection.
!>
!> This module implements the quadratic objective \(\phi\) used by the
!> projector to find the closest point on the DROP level-set surface.
!> The level-set constraint is owned and evaluated by the projector.
module moist_cavity_drop_objective_phi
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   implicit none

   integer, parameter :: ndim = 3

   !> Quadratic closest-point projection objective phi
   !>
   !> The objective is
   !> $$
   !> \phi(\mathbf{r}) = \frac{w_a}{2}\|\mathbf{r} - \mathbf{r}^\circ\|^2
   !> $$
   !> where $\mathbf{r}^\circ$ is the anchor point and $w_a$ is `param%phi_alpha`.
   type :: moist_cavity_drop_objective_phi_type
      !> Parameters for the DROP cavity
      type(moist_cavity_drop_parameters_type) :: param

      !> Number of molecular centers, used for nuclear derivative dimensions
      integer :: ncenters = 0

   contains
      !> Setup
      procedure :: set_parameters => phi_type_set_parameters
      procedure :: set_input => phi_type_set_input
      !> Value.
      procedure :: f0 => phi0
      !> Point derivatives
      procedure :: f1_r => phi1_r
      procedure :: f2_rr => phi2_rr
      procedure :: f3_rrr => phi3_rrr
      procedure :: f4_rrrr => phi4_rrrr
      !> Combined value and derivatives
      procedure :: f012_r => phi012_r
      !> Nuclear derivatives
      procedure :: f1_rA => phi1_rA
      procedure :: f2_rArB => phi2_rArB
      !> Mixed derivatives
      procedure :: f2_r_rA => phi2_r_rA
   end type moist_cavity_drop_objective_phi_type

   public :: moist_cavity_drop_objective_phi_type

contains

   !> Set objective parameters
   !> @param[inout] self Objective instance
   !> @param[in]    param DROP cavity parameters
   subroutine phi_type_set_parameters(self, param)
      class(moist_cavity_drop_objective_phi_type), intent(inout) :: self
      type(moist_cavity_drop_parameters_type), intent(in) :: param

      self%param = param
   end subroutine phi_type_set_parameters

   !> Record molecular input dimensions
   !> @param[inout] self Objective instance
   !> @param[in]    mol Molecular structure
   !> @param[in]    radii Atomic radii
   subroutine phi_type_set_input(self, mol, radii)
      class(moist_cavity_drop_objective_phi_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)

      self%ncenters = size(radii)
   end subroutine phi_type_set_input

   !> Compute phi
   !> @param[inout] self Objective instance
   !> @param[in]    pt Evaluation point
   !> @param[in]    anch Anchor coordinates
   !> @param[in]    owner Anchor owner index
   !> @returns      val phi value
   function phi0(self, pt, anch, owner) result(val)
      class(moist_cavity_drop_objective_phi_type), intent(inout) :: self
      real(wp), intent(in) :: pt(ndim)
      real(wp), intent(in) :: anch(ndim)
      integer, intent(in) :: owner
      real(wp) :: val
      real(wp) :: diff(ndim)

      diff = pt - anch
      val = 0.5_wp*self%param%phi_alpha*sum(diff*diff)
   end function phi0

   !> Compute the spatial gradient of phi
   !> @param[inout] self Objective instance
   !> @param[in]    pt Evaluation point
   !> @param[in]    anch Anchor coordinates
   !> @param[in]    owner Anchor owner index
   !> @returns      gradient Spatial gradient
   function phi1_r(self, pt, anch, owner) result(gradient)
      class(moist_cavity_drop_objective_phi_type), intent(inout) :: self
      real(wp), intent(in) :: pt(ndim)
      real(wp), intent(in) :: anch(ndim)
      integer, intent(in) :: owner
      real(wp) :: gradient(ndim)

      gradient = self%param%phi_alpha*(pt - anch)
   end function phi1_r

   !> Compute the spatial Hessian of phi
   !> @param[inout] self Objective instance
   !> @param[in]    pt Evaluation point
   !> @param[in]    anch Anchor coordinates
   !> @param[in]    owner Anchor owner index
   !> @returns      hessian Spatial Hessian
   function phi2_rr(self, pt, anch, owner) result(hessian)
      class(moist_cavity_drop_objective_phi_type), intent(inout) :: self
      real(wp), intent(in) :: pt(ndim)
      real(wp), intent(in) :: anch(ndim)
      integer, intent(in) :: owner
      real(wp) :: hessian(ndim, ndim)
      integer :: i

      hessian = 0.0_wp
      do i = 1, ndim
         hessian(i, i) = self%param%phi_alpha
      end do
   end function phi2_rr

   !> Compute phi, spatial gradient, and spatial Hessian
   !> @param[inout] self Objective instance
   !> @param[in]    pt Evaluation point
   !> @param[in]    anch Anchor coordinates
   !> @param[in]    owner Anchor owner index
   !> @param[out]   val phi value
   !> @param[out]   gradient Optional spatial gradient
   !> @param[out]   hessian Optional spatial Hessian
   subroutine phi012_r(self, pt, anch, owner, val, gradient, hessian)
      class(moist_cavity_drop_objective_phi_type), intent(inout) :: self
      real(wp), intent(in) :: pt(ndim)
      real(wp), intent(in) :: anch(ndim)
      integer, intent(in) :: owner
      real(wp), intent(out) :: val
      real(wp), intent(out), optional :: gradient(ndim)
      real(wp), intent(out), optional :: hessian(ndim, ndim)
      real(wp) :: diff(ndim)
      integer :: i

      diff = pt - anch
      val = 0.5_wp*self%param%phi_alpha*sum(diff*diff)

      if (present(gradient)) gradient = self%param%phi_alpha*diff

      if (present(hessian)) then
         hessian = 0.0_wp
         do i = 1, ndim
            hessian(i, i) = self%param%phi_alpha
         end do
      end if
   end subroutine phi012_r

   !> Compute the gradient of phi with respect to nuclear positions
   !> @param[in] self Objective instance
   !> @param[in] pt Evaluation point
   !> @param[in] anch Anchor coordinates
   !> @param[in] ownid Owner atom index
   !> @returns   gradient Nuclear gradient
   pure function phi1_rA(self, pt, anch, ownid) result(gradient)
      class(moist_cavity_drop_objective_phi_type), intent(in) :: self
      real(wp), intent(in) :: pt(ndim)
      real(wp), intent(in) :: anch(ndim)
      integer, intent(in) :: ownid
      real(wp), allocatable :: gradient(:, :)

      allocate (gradient(ndim, self%ncenters), source=0.0_wp)
      if (ownid < 1 .or. ownid > self%ncenters) return

      gradient(:, ownid) = -self%param%phi_alpha*(pt - anch)
   end function phi1_rA

   !> Compute the Hessian of phi with respect to nuclear positions
   !> @param[in] self Objective instance
   !> @param[in] pt Evaluation point
   !> @param[in] anch Anchor coordinates
   !> @param[in] ownid Owner atom index
   !> @returns   hessian Nuclear Hessian
   pure function phi2_rArB(self, pt, anch, ownid) result(hessian)
      class(moist_cavity_drop_objective_phi_type), intent(in) :: self
      real(wp), intent(in) :: pt(ndim)
      real(wp), intent(in) :: anch(ndim)
      integer, intent(in) :: ownid
      real(wp), allocatable :: hessian(:, :, :, :)
      integer :: axis

      allocate (hessian(ndim, ndim, self%ncenters, self%ncenters), source=0.0_wp)
      if (ownid < 1 .or. ownid > self%ncenters) return

      do axis = 1, ndim
         hessian(axis, axis, ownid, ownid) = self%param%phi_alpha
      end do
   end function phi2_rArB

   !> Compute the mixed spatial-nuclear Hessian of phi
   !> @param[in] self Objective instance
   !> @param[in] pt Evaluation point
   !> @param[in] anch Anchor coordinates
   !> @param[in] ownid Owner atom index
   !> @returns   deriv Mixed derivative
   pure function phi2_r_rA(self, pt, anch, ownid) result(deriv)
      class(moist_cavity_drop_objective_phi_type), intent(in) :: self
      real(wp), intent(in) :: pt(ndim)
      real(wp), intent(in) :: anch(ndim)
      integer, intent(in) :: ownid
      real(wp), allocatable :: deriv(:, :, :)
      integer :: axis

      allocate (deriv(ndim, ndim, self%ncenters), source=0.0_wp)
      if (ownid < 1 .or. ownid > self%ncenters) return

      do axis = 1, ndim
         deriv(axis, axis, ownid) = -self%param%phi_alpha
      end do
   end function phi2_r_rA

   !> Compute the third spatial derivative of phi
   !> @param[in] self Objective instance
   !> @param[in] pt Evaluation point
   !> @param[in] anch Anchor coordinates
   !> @param[in] owner Anchor owner index
   !> @returns   deriv Zero third derivative tensor
   pure function phi3_rrr(self, pt, anch, owner) result(deriv)
      class(moist_cavity_drop_objective_phi_type), intent(in) :: self
      real(wp), intent(in) :: pt(ndim)
      real(wp), intent(in) :: anch(ndim)
      integer, intent(in) :: owner
      real(wp) :: deriv(ndim, ndim, ndim)

      deriv = 0.0_wp
   end function phi3_rrr

   !> Compute the fourth spatial derivative of phi
   !> @param[in] self Objective instance
   !> @param[in] pt Evaluation point
   !> @param[in] anch Anchor coordinates
   !> @param[in] owner Anchor owner index
   !> @returns   deriv Zero fourth derivative tensor
   pure function phi4_rrrr(self, pt, anch, owner) result(deriv)
      class(moist_cavity_drop_objective_phi_type), intent(in) :: self
      real(wp), intent(in) :: pt(ndim)
      real(wp), intent(in) :: anch(ndim)
      integer, intent(in) :: owner
      real(wp) :: deriv(ndim, ndim, ndim, ndim)

      deriv = 0.0_wp
   end function phi4_rrrr

end module moist_cavity_drop_objective_phi
