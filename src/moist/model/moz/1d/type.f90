!> One-dimensional MOZ model state
module moist_model_moz_1d
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_model_type, only: solvation_model_type
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   implicit none(type, external)
   private
   public :: model_moz_1d_type, new_moz_1d_model

   !> Direct 1D MOZ model; radial theory is pending
   type, extends(solvation_model_type) :: model_moz_1d_type
      private
      !> Number of solute sites in the latest structure
      integer :: natom = 0
   contains
      procedure :: update => moz_1d_update
      procedure :: get_energy => moz_1d_energy
      procedure :: get_response => moz_1d_response
      procedure :: get_gradient => moz_1d_gradient
      procedure :: atom_count => moz_1d_atom_count
   end type model_moz_1d_type
contains
   !> Attach a borrowed run context
   !>
   !> @param[out] self Model
   !> @param[in] ctx Run context
   subroutine new_moz_1d_model(self, ctx)
      !> Self
      type(model_moz_1d_type), intent(out) :: self
      !> Ctx
      type(moist_context_type), intent(in), target :: ctx
      self%ctx => ctx
   end subroutine new_moz_1d_model

   !> Record the solute site count
   !>
   !> @param[in,out] self Model
   !> @param[in] mol Solute structure
   !> @param[out] error Invalid structure
   subroutine moz_1d_update(self, mol, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Mol
      class(structure_type), intent(in) :: mol
      !> Error
      type(error_type), allocatable, intent(out) :: error
      self%natom = 0
      if (.not. associated(self%ctx)) then
         call fatal_error(error, "Construct the 1D MOZ model before update")
         return
      end if
      if (mol%nat < 1) then
         call fatal_error(error, "1D MOZ requires at least one solute atom")
         return
      end if
      self%natom = mol%nat
   end subroutine moz_1d_update

   !> Number of sites in the latest valid structure
   !>
   !> @param[in] self Model
   function moz_1d_atom_count(self) result(natom)
      !> Self
      class(model_moz_1d_type), intent(in) :: self
      !> Natom
      integer :: natom
      natom = self%natom
   end function moz_1d_atom_count

   !> Report pending 1D MOZ energy theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] energy Energy
   !> @param[out] error Error
   subroutine moz_1d_energy(self, coupling, energy, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Energy
      real(wp), intent(inout) :: energy
      !> Error
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "1D MOZ energy is not implemented")
   end subroutine moz_1d_energy

   !> Report pending 1D MOZ response theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] response Response
   !> @param[out] error Error
   subroutine moz_1d_response(self, coupling, response, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Response
      type(response_type), intent(inout) :: response
      !> Error
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "1D MOZ response is not implemented")
   end subroutine moz_1d_response

   !> Report pending 1D MOZ gradient theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] response Response
   !> @param[in,out] gradient Gradient
   !> @param[out] error Error
   subroutine moz_1d_gradient(self, coupling, response, gradient, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Response
      type(response_type), intent(inout) :: response
      !> Gradient
      real(wp), intent(inout) :: gradient(:, :)
      !> Error
      type(error_type), allocatable, intent(out) :: error
      call fatal_error(error, "1D MOZ gradient is not implemented")
   end subroutine moz_1d_gradient
end module moist_model_moz_1d
