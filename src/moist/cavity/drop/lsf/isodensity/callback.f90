!> Callback-backed isodensity level set function for DROP
!>
!> This experimental LSF delegates value, spatial gradient, spatial Hessian,
!> and third spatial derivative evaluation to a C callback
module moist_cavity_drop_lsf_isodensity_callback
   use, intrinsic :: iso_c_binding, only: c_double, c_int, c_funptr, c_null_funptr, c_ptr, c_null_ptr, &
                            c_associated, c_f_procpointer, c_loc
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_base_update, lsf_candidate_space_user
   use moist_output_format, only: format_string
   implicit none (type, external)
   private

   integer, parameter :: ndim = 3

   public :: moist_cavity_drop_lsf_isodensity_callback_type
   public :: isodensity_lsf_callback

   abstract interface
      !> C callback for LSF value, gradient, and optional higher derivatives
      !>
      !> Value and gradient are always requested.  ``hess`` and ``third`` are
      !> passed as raw pointers that are NULL when that order is not required
      !> (driven by the cavity's ``set_max_deriv``): a NULL pointer signals the
      !> callee to skip computing -- not just writing -- that derivative, so the
      !> expensive density Hessian/third derivative is never evaluated during the
      !> value+gradient-only projection phase.
      !>
      !> @param[in]  context  User-owned callback context
      !> @param[in]  point    Evaluation point in Bohr
      !> @param[out] value    LSF value
      !> @param[out] grad     Spatial gradient dS/dr
      !> @param[out] hess     Spatial Hessian d2S/drdr   (double[3][3] or NULL)
      !> @param[out] third    Spatial third deriv         (double[3][3][3] or NULL)
      !> @returns             0 on success, nonzero host status on failure
      function isodensity_lsf_callback(context, point, value, grad, hess, third) result(status) bind(C)
         import :: c_double, c_int, c_ptr
         implicit none (type, external)
         type(c_ptr), value :: context
         real(c_double), intent(in) :: point(3)
         real(c_double), intent(out) :: value
         real(c_double), intent(out) :: grad(3)
         type(c_ptr), value :: hess
         type(c_ptr), value :: third
         integer(c_int) :: status
      end function isodensity_lsf_callback
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
      procedure, public :: f0 => lsf_f0
      procedure, public :: f012_r => lsf_f012_r
      procedure, public :: f3_rrr => lsf_f3_rrr
      procedure, public :: f3_rr_rA => lsf_f3_rr_rA
      procedure, public :: screening_offset => lsf_screening_offset
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

      !> The callback is globally evaluable and holds no per-atom data, so
      !> candidate ids are never translated.
      self%candidate_space = lsf_candidate_space_user

      self%radius_dependent = .false.

      self%callback_ptr = callback_ptr
      self%context = context
      if (present(scale)) self%scale = scale
   end subroutine lsf_new

   !> Bind molecular geometry for the inherited base state
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    mol   Molecular structure
   !> @param[in]    radii Per-atom radii
   subroutine lsf_update(self, mol, radii)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)

      call lsf_base_update(self, mol, radii)
   end subroutine lsf_update

   !> Evaluate and cache callback data at one point.
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    point Evaluation point in Bohr
   !> @param[out]   error Host callback failure at this point
   subroutine lsf_prepare(self, point, error)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      type(error_type), allocatable, intent(out) :: error

      procedure(isodensity_lsf_callback), pointer :: callback
      real(c_double) :: c_point(3), c_value, c_grad(3)
      real(c_double), target :: c_hess(3, 3), c_third(3, 3, 3)
      type(c_ptr) :: p_hess, p_third
      logical :: want_hess, want_third
      integer(c_int) :: status

      ! Only request the (expensive) density Hessian/third derivative for the
      ! orders the cavity actually needs.  The projection's value+gradient phase
      ! sets max_deriv=1, so a NULL hess/third pointer tells the callback to skip
      ! computing them entirely rather than evaluating and discarding them.
      want_hess = self%max_deriv >= 2
      want_third = self%max_deriv >= 3
      p_hess = c_null_ptr
      p_third = c_null_ptr
      if (want_hess) p_hess = c_loc(c_hess)
      if (want_third) p_third = c_loc(c_third)

      call c_f_procpointer(self%callback_ptr, callback)
      self%point = point
      !> Value and gradient are mandatory in the callback ABI; the higher orders
      !> are exactly those whose pointer was non-NULL.  Recorded so an accessor
      !> asked for more aborts instead of returning the zeros below.
      self%prepared_deriv = 1
      if (want_hess) self%prepared_deriv = 2
      if (want_third) self%prepared_deriv = 3
      c_point = real(point, c_double)
      status = callback(self%context, c_point, c_value, c_grad, p_hess, p_third)
      if (status /= 0_c_int) then
         ! Substitute state; numbers are not a result; the caller aborts on `error`
         self%value = 1.0_wp
         self%grad = 0.0_wp
         self%grad(1) = 1.0_wp
         self%hess = 0.0_wp
         self%third = 0.0_wp
         ! The substitute state is not a derivative jet of anything
         self%prepared_deriv = -1
         ! Carry the point in the diagnostic
         call fatal_error(error, "External LSF evaluation failed with status "// &
                          format_string(int(status), "(i0)")//" at point ("// &
                          trim(adjustl(format_string(point(1), "(es13.6)")))//" "// &
                          trim(adjustl(format_string(point(2), "(es13.6)")))//" "// &
                          trim(adjustl(format_string(point(3), "(es13.6)")))// &
                          " ) Bohr. The cavity build was aborted; no cavity data are valid.")
         return
      end if
      self%value = self%scale*real(c_value, wp)
      self%grad = self%scale*real(c_grad, wp)
      if (want_hess) then
         self%hess = self%scale*real(c_hess, wp)
      else
         self%hess = 0.0_wp
      end if
      if (want_third) then
         self%third = self%scale*real(c_third, wp)
      else
         self%third = 0.0_wp
      end if
   end subroutine lsf_prepare

   !> Evaluate callback data; candidate lists are ignored for true density LSFs
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point in Bohr
   !> @param[in]    candidate_indices Ignored atom candidates
   !> @param[out]   error             Host callback failure at this point
   subroutine lsf_prepare_subset(self, point, candidate_indices, error)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(inout) :: self
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: candidate_indices(:)
      type(error_type), allocatable, intent(out) :: error

      if (size(candidate_indices) < 0) return
      call self%prepare(point, error)
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

   !> Number of active atoms. True-density callbacks are not atom screened
   pure function lsf_active_count(self) result(n)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      integer :: n

      if (self%ncenters < 0) then
         n = -1
         return
      end if
      n = 0
   end function lsf_active_count

   !> Active atom lookup. Undefined for zero active atoms, returns zero sentinel
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
   subroutine lsf_f0(self, val)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      real(wp), intent(out) :: val

      val = self%value
   end subroutine lsf_f0

   !> Return cached LSF value, gradient, and Hessian.
   !>
   !> @param[in]  self    LSF instance
   !> @param[out] lsf0    Optional LSF value
   !> @param[out] lsf1_r  Optional spatial gradient
   !> @param[out] lsf2_rr Optional spatial Hessian
   subroutine lsf_f012_r(self, lsf0, lsf1_r, lsf2_rr)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)

      if (present(lsf0)) lsf0 = self%value
      if (present(lsf1_r)) lsf1_r(:) = self%grad(:)
      if (present(lsf2_rr)) then
         call self%require_deriv(2, "f012_r(lsf2_rr)")
         lsf2_rr(:, :) = self%hess(:, :)
      end if
   end subroutine lsf_f012_r

   !> Return cached lower derivatives and third spatial derivative
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf0      Optional LSF value
   !> @param[out] lsf1_r    Optional spatial gradient
   !> @param[out] lsf2_rr   Optional spatial Hessian
   !> @param[out] lsf3_rrr  Spatial third derivative
   subroutine lsf_f3_rrr(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf0
      real(wp), intent(out), optional :: lsf1_r(:)
      real(wp), intent(out), optional :: lsf2_rr(:, :)
      real(wp), intent(out) :: lsf3_rrr(:, :, :)

      call self%require_deriv(3, "f3_rrr")
      call self%f012_r(lsf0, lsf1_r, lsf2_rr)
      lsf3_rrr(:, :, :) = self%third
   end subroutine lsf_f3_rrr

   !> Return zero nuclear derivative placeholders
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf1_rA     Optional nuclear gradient placeholder
   !> @param[out] lsf2_r_rA   Optional mixed second derivative placeholder
   !> @param[out] lsf3_rr_rA  Mixed third derivative placeholder
   subroutine lsf_f3_rr_rA(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      real(wp), intent(out) :: lsf3_rr_rA(:, :, :, :)

      if (present(lsf1_rA)) lsf1_rA(:, :) = 0.0_wp
      if (present(lsf2_r_rA)) lsf2_r_rA(:, :, :) = 0.0_wp
      lsf3_rr_rA(:, :, :, :) = 0.0_wp
   end subroutine lsf_f3_rr_rA

   !> Density callback is globally evaluable; no atom-specific reach is needed
   !>
   !> Nothing about the callback decays with distance from an atom, so the offset
   !> is zero: the cavity cell grid gets no extra shell and the candidate lists
   !> this LSF is handed are ignored anyway.
   !>
   !> @param[in] self   LSF instance
   !> @param[in] radius Atom radius (Bohr)
   !> @returns          Radial offset from the atom surface (Bohr), always zero
   pure function lsf_screening_offset(self, radius) result(d)
      class(moist_cavity_drop_lsf_isodensity_callback_type), intent(in) :: self
      !> Atom radius (Bohr)
      real(wp), intent(in) :: radius
      !> Radial offset from the atom surface (Bohr)
      real(wp) :: d

      if (self%ncenters < 0) then
         d = radius
      else
         d = 0.0_wp*radius
      end if
   end function lsf_screening_offset

end module moist_cavity_drop_lsf_isodensity_callback
