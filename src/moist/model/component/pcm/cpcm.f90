!> CPCM (Conductor-like Polarizable Continuum Model) implementation
!> This module provides the CPCM variant of PCM with its specific dielectric
!> scaling (f epsilon = ( epsilon -1)/ epsilon ). Matrix assembly is delegated to the cavity type.
module moist_model_component_pcm_cpcm
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_type, only: cavity_type, wavefunction_type
   use moist_model_component_pcm_type, only: pcm_base, solver_type, &
      & potential_source
   implicit none
   private

   public :: cpcm
   public :: new_cpcm

   !> CPCM (Conductor-like Polarizable Continuum Model) variant
   !> Uses f epsilon = ( epsilon -1)/ epsilon scaling. Matrix assembly is delegated to the cavity type.
   type, extends(pcm_base) :: cpcm
   end type cpcm

contains

   !> Constructor for CPCM variant
   !> Sets f epsilon = ( epsilon -1)/ epsilon and configures solver and potential source.
   subroutine new_cpcm(self, epsilon, solver, phi_source, external_matrix, error)
      !> CPCM instance to initialize
      type(cpcm), intent(out) :: self
      !> Dielectric constant
      real(wp), intent(in) :: epsilon
      !> Optional: linear solver type
      integer, intent(in), optional :: solver
      !> Optional: potential source strategy
      integer, intent(in), optional :: phi_source
      !> Optional: external pre-computed matrix (ngrid, ngrid)
      real(wp), intent(in), optional :: external_matrix(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      ! Set dielectric properties
      self%epsilon = epsilon
      self%feps = (epsilon - 1.0_wp)/epsilon  ! CPCM formula

      ! Set solver type
      if (present(solver)) then
         self%solver = solver
      else
         self%solver = solver_type%cholesky
      end if

      ! Set potential source
      if (present(phi_source)) then
         self%phi_source = phi_source
      else
         self%phi_source = potential_source%charges
      end if

      ! Handle external matrix
      if (present(external_matrix)) then
         call self%set_external_matrix(external_matrix)
      end if

      ! Set component name
      self%name = "CPCM"

   end subroutine new_cpcm

end module moist_model_component_pcm_cpcm
