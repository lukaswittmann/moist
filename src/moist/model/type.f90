!> Domain-independent solvation model interface
module moist_model_type
   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_channels_response, only: response_type
   use moist_channels_coupling, only: coupling_type

   implicit none(type, external)
   private

   public :: solvation_model_type

   !> Abstract base solvation model
   type, abstract :: solvation_model_type
      !> Borrowed run context (verbosity/debug/timer); set at construction,
      !> owned by the top-level caller, never allocated or freed by the model
      type(moist_context_type), pointer :: ctx => null()

   contains

      procedure(update_model), deferred :: update
      procedure(get_model_energy), deferred :: get_energy
      procedure(get_model_response), deferred :: get_response
      procedure(get_model_gradient), deferred :: get_gradient

   end type solvation_model_type

   abstract interface

      !> Update the solvation model with the current molecular structure
      !>
      !> - calculates all structure-dependent properties
      subroutine update_model(self, mol, error)
         import solvation_model_type, structure_type, error_type
         implicit none(type, external)
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Molecular structure data
         class(structure_type), intent(in) :: mol
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_model

      !> Evaluate the solvation energy
      subroutine get_model_energy(self, coupling, energy, error)
         import solvation_model_type, structure_type, wp, error_type, coupling_type
         implicit none(type, external)
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(inout), target :: coupling
         !> Solvation energy
         real(wp), intent(inout) :: energy
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_energy

      !> Get the solvation response (only for self-consistent models)
      subroutine get_model_response(self, coupling, response, error)
         import solvation_model_type, structure_type, wp, error_type, response_type, coupling_type
         implicit none(type, external)
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(inout), target :: coupling
         !> Solvation response for the component
         type(response_type), intent(inout) :: response
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_response

      !> Get the solvation energy gradient and the host part of the phase
      !>
      !> `response` is cleared on entry and returns what the host contracts
      !> with its own geometry derivatives (potential adjoint, Gaussian amplitudes)
      subroutine get_model_gradient(self, coupling, response, gradient, error)
         import solvation_model_type, structure_type, wp, error_type, response_type, coupling_type
         implicit none(type, external)
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(inout), target :: coupling
         !> Host part of the gradient phase
         type(response_type), intent(inout) :: response
         !> Solvation gradient
         real(wp), intent(inout) :: gradient(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_gradient

   end interface

end module moist_model_type
