!> Data envelopes exchanged between a host program and moist
!>
!> The directions treat a missing array differently:
!>
!> * A missing input is a host mistake; consumers assert what they need at
!>   the point of use with [[require_channel]]
!> * A missing output reflects correct physics; a cavity with, e.g.,
!>   field-independent geometry has no level set response; or a model without
!>   GOSTSHYP has no Gaussian amplitudes
module moist_channels
   use mctc_env, only: wp, error_type, fatal_error

   implicit none
   private

   public :: coupling_type, response_type
   public :: coupling_channel_type, response_channel_type
   public :: electrostatic_coupling_type, gostshyp_coupling_type
   public :: electrostatic_response_type, lsf_response_type, gostshyp_response_type
   public :: require_channel

   !> Assert that a host-supplied coupling channel is present and correctly shaped
   !>
   !> Consumers call this where they read a channel. An absent input is a host
   !> mistake and must be reported, never silently read as zero.
   interface require_channel
      module procedure :: require_channel_vector
      module procedure :: require_channel_matrix
      module procedure :: require_channel_tensor3
   end interface require_channel

   !* ============================================================================== *!
   !*                              Channel base types                                *!
   !* ============================================================================== *!

   !> Base type for one group of host-supplied coupling channels
   !>
   !> The binding fixes a common interface across the groups; the bodies cannot
   !> be shared, because Fortran offers no way to enumerate a type's components.
   type, abstract :: coupling_channel_type
   contains
      !> Deallocate every channel in this group
      procedure(clear_coupling_channel), deferred :: clear
   end type coupling_channel_type

   !> Base type for one group of channels returned to the host
   type, abstract :: response_channel_type
   contains
      !> Deallocate every channel in this group
      procedure(clear_response_channel), deferred :: clear
   end type response_channel_type

   abstract interface

      !> Deallocate every channel in a coupling group
      subroutine clear_coupling_channel(self)
         import :: coupling_channel_type
         !> Coupling group to clear
         class(coupling_channel_type), intent(inout) :: self
      end subroutine clear_coupling_channel

      !> Deallocate every channel in a response group
      subroutine clear_response_channel(self)
         import :: response_channel_type
         !> Response group to clear
         class(response_channel_type), intent(inout) :: self
      end subroutine clear_response_channel

   end interface

   !* ============================================================================== *!
   !*                        Host -> MOIST coupling channels                         *!
   !* ============================================================================== *!

   !> Electrostatic channels supplied by the host for one coupling step
   type, extends(coupling_channel_type) :: electrostatic_coupling_type
      !> Number of electrons for each atom (nat, spin)
      real(wp), allocatable :: qat(:, :)
      !> Molecular potential trace on the electrostatic cavity (ngrid)
      real(wp), allocatable :: phi(:)
      !> Direct host trace-geometry weights for Gaussian widths (ngrid)
      real(wp), allocatable :: w_xi(:)
      !> Direct host trace-geometry weights for switch factors (ngrid)
      real(wp), allocatable :: w_f(:)
      !> Direct host trace-geometry weights for positions (3, ngrid)
      real(wp), allocatable :: w_xyz(:, :)
      !> Direct host trace-geometry weights for normals (3, ngrid)
      real(wp), allocatable :: w_normal(:, :)
      !> Charge-weighted electronic field on the electrostatic cavity (3, ngrid)
      !  (distinct from `w_xyz`, which is an adjoint with respect to surface positions)
      real(wp), allocatable :: qefield(:, :)
   contains
      !> Clear all supplied electrostatic channels
      procedure :: clear => clear_electrostatic_coupling
   end type electrostatic_coupling_type

   !> GOSTSHYP channels supplied by the host for one coupling step
   type, extends(coupling_channel_type) :: gostshyp_coupling_type
      !> Moments of the solute density against an unnormalized Gaussian
      !> `G_i = exp(-w_i |r - r_i|^2)` centered on a cavity grid point
      !>
      !> The width `w_i` is the model's own, so the host must read it back
      !> from the cavity before forming these
      !>
      !> AO-basis three-center integrals contracted with the density, needs
      !> QM density and integrals
      !>
      !> `<G_i>` (ngrid)
      real(wp), allocatable :: gt(:)
      !> `<(r - r_i) G_i>` (3, ngrid)
      real(wp), allocatable :: pt(:, :)
      !> `<(r - r_i)(r - r_i) G_i>` (3, 3, ngrid) (symmetric in the first two)
      real(wp), allocatable :: mt(:, :, :)
      !> `<(r - r_i) |r - r_i|^2 G_i>` (3, ngrid)
      real(wp), allocatable :: rt(:, :)
   contains
      !> Clear all supplied GOSTSHYP moments
      procedure :: clear => clear_gostshyp_coupling
   end type gostshyp_coupling_type

   !> QM/solute data supplied to solvation models for one coupling step
   !>
   !> Every channel is optional at the type level; a consumer that needs one
   !> validates it at the point of use with [[require_channel]] rather than
   !> silently treating an absent array as zero.
   type :: coupling_type
      !> Electrostatic channels
      type(electrostatic_coupling_type) :: electrostatics
      !> GOSTSHYP channels
      type(gostshyp_coupling_type) :: gostshyp
   contains
      !> Clear all supplied coupling arrays
      procedure :: clear => clear_coupling
   end type coupling_type

   !* ============================================================================== *!
   !*                        MOIST -> host response channels                         *!
   !* ============================================================================== *!

   !> Electrostatic response returned to the host
   type, extends(response_channel_type) :: electrostatic_response_type
      !> Surface charge `q_i`, which equals the adjoint `dE/dphi_i` by
      !> stationarity. The host contracts it with its own potential integrals,
      !> `F_uv += sum_i q_i V_uv(r_i)`.
      !>
      !> This is an accumulator over components; only charge-like contributions
      !> may be added here, so that the name stays true of the sum.
      real(wp), allocatable :: surface_charge(:)
   contains
      !> Clear the electrostatic response
      procedure :: clear => clear_electrostatic_response
   end type electrostatic_response_type

   !> Cavity level set response returned to the host
   type, extends(response_channel_type) :: lsf_response_type
      !> Adjoint weights for cavity level set values (ngrid)
      real(wp), allocatable :: w_value(:)
      !> Adjoint weights for cavity level set gradients (3, ngrid)
      real(wp), allocatable :: w_gradient(:, :)
      !> Adjoint weights for cavity level set Hessians (3, 3, ngrid)
      real(wp), allocatable :: w_hessian(:, :, :)
   contains
      !> Clear the level set response
      procedure :: clear => clear_lsf_response
   end type lsf_response_type

   !> GOSTSHYP response returned to the host
   type, extends(response_channel_type) :: gostshyp_response_type
      !> Amplitudes conjugate to the host's Gaussian integral blocks (ngrid)
      !> (counterpart of the `gostshyp` coupling moments)
      !> Host builds Fock contribution as plain sum over grid points,
      !>
      !>    F_uv += sum_i [ w_overlap(i) g_uv,i + w_normal_deriv(i) f_uv,i ]
      !>
      !> with `g_uv,i = <u|G_i|v>` and `f_uv,i = n_i . grad_r g_uv,i`
      !>
      !> Both signs are folded in here so the host never has to know the convention
      !>
      !> Grid points the model has switched off carry exactly zero, so the mask
      !> propagates without the host repeating it
      real(wp), allocatable :: w_overlap(:)
      !> See `w_overlap` (ngrid)
      real(wp), allocatable :: w_normal_deriv(:)
   contains
      !> Clear the GOSTSHYP response
      procedure :: clear => clear_gostshyp_response
   end type gostshyp_response_type

   !> Solvation response handed back to the host for one coupling step
   !>
   !> A channel left unallocated means the model has no contribution to it,
   !> which is a statement about the physics rather than an error: a cavity
   !> with field-independent geometry has no level set response, and a model
   !> without a GOSTSHYP component has no Gaussian amplitudes. Hosts that
   !> *require* a channel assert that at the API boundary.
   type :: response_type
      !> Electrostatic response
      type(electrostatic_response_type) :: electrostatics
      !> Cavity level set response
      type(lsf_response_type) :: lsf
      !> GOSTSHYP response
      type(gostshyp_response_type) :: gostshyp
   contains
      !> Clear all accumulated response arrays
      procedure :: clear => clear_response
   end type response_type

contains

   !> Clear all arrays supplied for one QM-solvation coupling step.
   subroutine clear_coupling(self)
      !> Coupling data to clear
      class(coupling_type), intent(inout) :: self

      call self%electrostatics%clear()
      call self%gostshyp%clear()

   end subroutine clear_coupling

   !> Clear the host-supplied electrostatic channels
   subroutine clear_electrostatic_coupling(self)
      !> Electrostatic coupling data to clear
      class(electrostatic_coupling_type), intent(inout) :: self

      if (allocated(self%qat)) deallocate (self%qat)
      if (allocated(self%phi)) deallocate (self%phi)
      if (allocated(self%w_xi)) deallocate (self%w_xi)
      if (allocated(self%w_f)) deallocate (self%w_f)
      if (allocated(self%w_xyz)) deallocate (self%w_xyz)
      if (allocated(self%w_normal)) deallocate (self%w_normal)
      if (allocated(self%qefield)) deallocate (self%qefield)

   end subroutine clear_electrostatic_coupling

   !> Clear the host-supplied GOSTSHYP moments
   subroutine clear_gostshyp_coupling(self)
      !> GOSTSHYP coupling data to clear
      class(gostshyp_coupling_type), intent(inout) :: self

      if (allocated(self%gt)) deallocate (self%gt)
      if (allocated(self%pt)) deallocate (self%pt)
      if (allocated(self%mt)) deallocate (self%mt)
      if (allocated(self%rt)) deallocate (self%rt)

   end subroutine clear_gostshyp_coupling

   !> Clear all accumulated response arrays
   subroutine clear_response(self)
      !> Response data to clear
      class(response_type), intent(inout) :: self

      call self%electrostatics%clear()
      call self%lsf%clear()
      call self%gostshyp%clear()

   end subroutine clear_response

   !> Clear the electrostatic response
   subroutine clear_electrostatic_response(self)
      !> Electrostatic response to clear
      class(electrostatic_response_type), intent(inout) :: self

      if (allocated(self%surface_charge)) deallocate (self%surface_charge)

   end subroutine clear_electrostatic_response

   !> Clear the level set response
   subroutine clear_lsf_response(self)
      !> level set response to clear
      class(lsf_response_type), intent(inout) :: self

      if (allocated(self%w_value)) deallocate (self%w_value)
      if (allocated(self%w_gradient)) deallocate (self%w_gradient)
      if (allocated(self%w_hessian)) deallocate (self%w_hessian)

   end subroutine clear_lsf_response

   !> Clear the GOSTSHYP response
   subroutine clear_gostshyp_response(self)
      !> GOSTSHYP response to clear
      class(gostshyp_response_type), intent(inout) :: self

      if (allocated(self%w_overlap)) deallocate (self%w_overlap)
      if (allocated(self%w_normal_deriv)) deallocate (self%w_normal_deriv)

   end subroutine clear_gostshyp_response

   !> Require a host-supplied vector channel of exactly `n` elements
   !>
   !> @param[in]  array Channel to validate
   !> @param[in]  name  Channel name used in diagnostics
   !> @param[in]  n     Expected number of elements
   !> @param[out] error Error handling
   subroutine require_channel_vector(array, name, n, error)
      !> Channel to validate
      real(wp), allocatable, intent(in) :: array(:)
      !> Channel name used in diagnostics
      character(len=*), intent(in) :: name
      !> Expected number of elements
      integer, intent(in) :: n
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(array)) then
         call missing_channel(name, error)
         return
      end if
      if (size(array) /= n) then
         call misshaped_channel(name, [size(array)], [n], error)
      end if

   end subroutine require_channel_vector

   !> Require a host-supplied rank-2 channel of shape `(d1, d2)`
   !>
   !> @param[in]  array Channel to validate
   !> @param[in]  name  Channel name used in diagnostics
   !> @param[in]  d1    Expected extent of the first dimension
   !> @param[in]  d2    Expected extent of the second dimension
   !> @param[out] error Error handling
   subroutine require_channel_matrix(array, name, d1, d2, error)
      !> Channel to validate
      real(wp), allocatable, intent(in) :: array(:, :)
      !> Channel name used in diagnostics
      character(len=*), intent(in) :: name
      !> Expected extent of the first dimension
      integer, intent(in) :: d1
      !> Expected extent of the second dimension
      integer, intent(in) :: d2
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(array)) then
         call missing_channel(name, error)
         return
      end if
      if (size(array, 1) /= d1 .or. size(array, 2) /= d2) then
         call misshaped_channel(name, shape(array), [d1, d2], error)
      end if

   end subroutine require_channel_matrix

   !> Require a host-supplied rank-3 channel of shape `(d1, d2, d3)`
   !>
   !> @param[in]  array Channel to validate
   !> @param[in]  name  Channel name used in diagnostics
   !> @param[in]  d1    Expected extent of the first dimension
   !> @param[in]  d2    Expected extent of the second dimension
   !> @param[in]  d3    Expected extent of the third dimension
   !> @param[out] error Error handling
   subroutine require_channel_tensor3(array, name, d1, d2, d3, error)
      !> Channel to validate
      real(wp), allocatable, intent(in) :: array(:, :, :)
      !> Channel name used in diagnostics
      character(len=*), intent(in) :: name
      !> Expected extent of the first dimension
      integer, intent(in) :: d1
      !> Expected extent of the second dimension
      integer, intent(in) :: d2
      !> Expected extent of the third dimension
      integer, intent(in) :: d3
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(array)) then
         call missing_channel(name, error)
         return
      end if
      if (size(array, 1) /= d1 .or. size(array, 2) /= d2 .or. size(array, 3) /= d3) then
         call misshaped_channel(name, shape(array), [d1, d2, d3], error)
      end if

   end subroutine require_channel_tensor3

   !> Report a coupling channel the host did not supply
   subroutine missing_channel(name, error)
      !> Channel name used in diagnostics
      character(len=*), intent(in) :: name
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "Required coupling channel '"//name// &
                       & "' was not supplied by the host")

   end subroutine missing_channel

   !> Report a coupling channel supplied with the wrong shape
   subroutine misshaped_channel(name, actual, expected, error)
      !> Channel name used in diagnostics
      character(len=*), intent(in) :: name
      !> Shape the host supplied
      integer, intent(in) :: actual(:)
      !> Shape the consumer requires
      integer, intent(in) :: expected(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "Coupling channel '"//name//"' has shape "// &
                       & shape_string(actual)//", expected "//shape_string(expected))

   end subroutine misshaped_channel

   !> Render a shape as `(d1, d2, ...)` for diagnostics
   function shape_string(extents) result(string)
      !> Extents to render
      integer, intent(in) :: extents(:)
      !> Rendered shape
      character(len=:), allocatable :: string

      character(len=32) :: buffer
      integer :: i

      string = "("
      do i = 1, size(extents)
         write (buffer, "(i0)") extents(i)
         if (i > 1) string = string//", "
         string = string//trim(buffer)
      end do
      string = string//")"

   end function shape_string

end module moist_channels
