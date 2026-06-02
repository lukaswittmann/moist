!> Callback-backed isodensity level-set function for DROP
!>
!> This experimental LSF delegates value, spatial gradient, spatial Hessian,
!> and third spatial derivative evaluation to a C callback.
module moist_cavity_drop_lsf_isodensity_callback
   use iso_c_binding, only: c_double, c_funptr, c_null_funptr, c_ptr, c_null_ptr, &
                            c_associated, c_f_procpointer
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   implicit none
   private

   integer, parameter :: ndim = 3

   public :: moist_cavity_drop_lsf_isodensity_callback_type
   public :: isodensity_lsf_callback

   abstract interface
      !> C callback for LSF value, gradient, Hessian, and third derivative at one point
      !>
      !> @param[in]  context  User-owned callback context
      !> @param[in]  point    Evaluation point in Bohr
      !> @param[out] value    LSF value
      !> @param[out] grad     Spatial gradient dS/dr
      !> @param[out] hess     Spatial Hessian d2S/drdr
      !> @param[out] third    Spatial third derivative d3S/drdrdr
      subroutine isodensity_lsf_callback(context, point, value, grad, hess, third) bind(C)
         import :: c_double, c_ptr
         type(c_ptr), value :: context
         real(c_double), intent(in) :: point(3)
         real(c_double), intent(out) :: value
         real(c_double), intent(out) :: grad(3)
         real(c_double), intent(out) :: hess(3, 3)
         real(c_double), intent(out) :: third(3, 3, 3)
      end subroutine isodensity_lsf_callback
   end interface

   !> Isodensity LSF implemented by a foreign callback.
   type, extends(moist_cavity_drop_lsf_type) :: moist_cavity_drop_lsf_isodensity_callback_type
      !> Raw C callback pointer
      type(c_funptr) :: callback_ptr = c_null_funptr
      !> User context passed through to the callback
      type(c_ptr) :: context = c_null_ptr
      !> Constant multiplier applied to value and spatial derivatives
      real(wp) :: scale = 100.0_wp
      !> Cached evaluation point in Bohr
      real(wp) :: point(ndim) = 0.0_wp
      !> Cached LSF value
      real(wp) :: value = 0.0_wp
      !> Cached spatial gradient
      real(wp) :: grad(ndim) = 0.0_wp
      !> Cached spatial Hessian
      real(wp) :: hess(ndim, ndim) = 0.0_wp
      !> Cached third spatial derivative
      real(wp) :: third(ndim, ndim, ndim) = 0.0_wp
      !> Highest requested derivative order
      integer :: max_deriv = 0
   contains
      procedure, public :: new => lsf_new
      procedure, public :: update => lsf_update
      procedure, public :: prepare => lsf_prepare
      procedure, public :: prepare_subset => lsf_prepare_subset
      procedure, public :: set_max_deriv => lsf_set_max_deriv
      procedure, public :: active_count => lsf_active_count
      procedure, public :: active_atom => lsf_active_atom
      procedure, public :: f0_screened => lsf_f0_screened
      procedure, public :: f012_r_screened => lsf_f012_r_screened
      procedure, public :: f3_rrr_screened => lsf_f3_rrr_screened
      procedure, public :: f3_rr_rA_screened => lsf_f3_rr_rA_screened
      procedure, public :: neighbor_cutoff => lsf_neighbor_cutoff
   end type moist_cavity_drop_lsf_isodensity_callback_type

contains

   !> Configure the callback pointer and context
   !>
   !> @param[inout] self         LSF instance
   !> @param[in]    callback_ptr C function pointer for density LSF evaluation
   !> @param[in]    context      User callback context
   !> @param[in]    scale        Constant LSF multiplier
   subroutine lsf_new(self, callback_ptr, context, scale)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: self
      type(c_funptr), intent(in) :: callback_ptr
      type(c_ptr), intent(in) :: context
      real(wp), intent(in), optional :: scale

      self%callback_ptr = callback_ptr
      self%context = context
      if (present(scale)) self%scale = scale
   end subroutine lsf_new

   !> Bind molecular geometry for the inherited base state.
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    mol   Molecular structure
   !> @param[in]    radii Per-atom radii
   subroutine lsf_update(self, mol, radii)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)

      self%mol = mol
      self%ncenters = mol%nat
      if (allocated(self%radii)) deallocate (self%radii)
      self%radii = radii
   end subroutine lsf_update

   !> Evaluate and cache callback data at one point.
   !> @param[inout] self  LSF instance
   !> @param[in]    point Evaluation point in Bohr
   subroutine lsf_prepare(self, point)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)

      procedure(isodensity_lsf_callback), pointer :: callback
      real(c_double) :: c_point(3), c_value, c_grad(3), c_hess(3, 3), c_third(3, 3, 3)

      ! TODO: Proper error propagration
      if (.not. c_associated(self%callback_ptr)) error stop "isodensity LSF callback is not associated"

      call c_f_procpointer(self%callback_ptr, callback)
      self%point = point
      c_point = real(point, c_double)
      call callback(self%context, c_point, c_value, c_grad, c_hess, c_third)
      self%value = self%scale*real(c_value, wp)
      self%grad = self%scale*real(c_grad, wp)
      self%hess = self%scale*real(c_hess, wp)
      self%third = self%scale*real(c_third, wp)
   end subroutine lsf_prepare

   !> Evaluate callback data; candidate lists are ignored for true density LSFs
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point in Bohr
   !> @param[in]    candidate_indices Ignored atom candidates
   subroutine lsf_prepare_subset(self, point, candidate_indices)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: candidate_indices(:)

      if (size(candidate_indices) < 0) return
      call self%prepare(point)
   end subroutine lsf_prepare_subset

   !> Record requested derivative order
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    n    Requested max derivative order
   subroutine lsf_set_max_deriv(self, n)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: self
      integer, intent(in) :: n

      self%max_deriv = max(0, n)
   end subroutine lsf_set_max_deriv

   !> Number of active atoms. True-density callbacks are not atom screened.
   pure function lsf_active_count(self) result(n)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      integer :: n

      if (self%ncenters < 0) then
         n = -1
         return
      end if
      n = 0
   end function lsf_active_count

   !> Active atom lookup. Undefined for zero active atoms, returns zero sentinel.
   pure function lsf_active_atom(self, i) result(idx)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      integer, intent(in) :: i
      integer :: idx

      if (self%ncenters < 0 .or. i < 0) then
         idx = -1
         return
      end if
      idx = 0
   end function lsf_active_atom

   !> Return cached LSF value
   !>
   !> @param[in]  self LSF instance
   !> @param[out] val  LSF value
   subroutine lsf_f0_screened(self, val)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      real(wp), intent(out) :: val

      val = self%value
   end subroutine lsf_f0_screened

   !> Return cached LSF value, gradient, and Hessian.
   !> @param[in]  self    LSF instance
   !> @param[out] lsf0    Optional LSF value
   !> @param[out] lsf1_r  Optional spatial gradient
   !> @param[out] lsf2_rr Optional spatial Hessian
   subroutine lsf_f012_r_screened(self, lsf0, lsf1_r, lsf2_rr)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)

      if (present(lsf0)) lsf0 = self%value
      if (present(lsf1_r)) lsf1_r(:) = self%grad(:)
      if (present(lsf2_rr)) lsf2_rr(:, :) = self%hess(:, :)
   end subroutine lsf_f012_r_screened

   !> Return cached lower derivatives and third spatial derivative.
   !> @param[in]  self      LSF instance
   !> @param[out] lsf0      Optional LSF value
   !> @param[out] lsf1_r    Optional spatial gradient
   !> @param[out] lsf2_rr   Optional spatial Hessian
   !> @param[out] lsf3_rrr  Spatial third derivative
   subroutine lsf_f3_rrr_screened(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)
      real(wp), allocatable, intent(out) :: lsf3_rrr(:, :, :)

      call self%f012_r_screened(lsf0, lsf1_r, lsf2_rr)
      allocate (lsf3_rrr(ndim, ndim, ndim), source=self%third)
   end subroutine lsf_f3_rrr_screened

   !> Return zero nuclear derivative placeholders
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf1_rA     Optional nuclear gradient placeholder
   !> @param[out] lsf2_r_rA   Optional mixed second derivative placeholder
   !> @param[out] lsf3_rr_rA  Mixed third derivative placeholder
   subroutine lsf_f3_rr_rA_screened(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      real(wp), allocatable, intent(out) :: lsf3_rr_rA(:, :, :, :)

      if (present(lsf1_rA)) lsf1_rA(:, :) = 0.0_wp
      if (present(lsf2_r_rA)) lsf2_r_rA(:, :, :) = 0.0_wp
      allocate (lsf3_rr_rA(ndim, ndim, ndim, self%ncenters), source=0.0_wp)
   end subroutine lsf_f3_rr_rA_screened

   !> Density callback is globally evaluable; no atom-specific reach is needed
   !>
   !> @param[in] self   LSF instance
   !> @param[in] radius Atom radius
   pure function lsf_neighbor_cutoff(self, radius) result(d)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      real(wp), intent(in) :: radius
      real(wp) :: d

      if (self%ncenters < 0) then
         d = radius
      else
         d = 0.0_wp*radius
      end if
   end function lsf_neighbor_cutoff

end module moist_cavity_drop_lsf_isodensity_callback
