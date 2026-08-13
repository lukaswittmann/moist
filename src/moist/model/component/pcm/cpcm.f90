!> CPCM (Conductor-like Polarizable Continuum Model) implementation
!> This module provides the CPCM variant of PCM with its specific dielectric
!> scaling (f epsilon = ( epsilon -1)/ epsilon )
module moist_model_component_pcm_cpcm
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_type, only: cavity_type
   use moist_channels, only: coupling_type
   use moist_model_component_pcm_type, only: solvation_model_component_pcm, solver_type, &
      & potential_source
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none (type, external)
   private

   public :: solvation_model_component_cpcm
   public :: new_component_cpcm

   !> CPCM (Conductor-like Polarizable Continuum Model) variant
   !> Uses f epsilon = ( epsilon -1)/ epsilon scaling
   type, extends(solvation_model_component_pcm) :: solvation_model_component_cpcm
   end type solvation_model_component_cpcm

contains

   !> Constructor for CPCM variant
   !> Sets f epsilon = ( epsilon -1)/ epsilon and configures solver and potential source.
   subroutine new_component_cpcm(self, ctx, epsilon, solver, phi_source, external_matrix, error)
      !> CPCM instance to initialize
      type(solvation_model_component_cpcm), intent(out) :: self
      !> Shared run context (verbosity/debug/timer); borrowed, must outlive self
      type(moist_context_type), intent(in), target :: ctx
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

      !> Borrow the shared run context (owns verbosity/debug/timer)
      self%ctx => ctx

      ! Set dielectric properties. Below eps = 1 the scaling factor turns
      ! negative (and diverges at eps = 0), so the model is undefined there.
      ! The `epsilon /= epsilon` test rejects a NaN input.
      if (epsilon < 1.0_wp .or. epsilon /= epsilon) then
         call fatal_error(error, &
            & "[new_component_cpcm] Dielectric constant must be >= 1")
         return
      end if
      self%epsilon = epsilon
      if (ieee_is_finite(epsilon)) then
         self%feps = (epsilon - 1.0_wp)/epsilon  ! CPCM formula
      else
         self%feps = 1.0_wp  ! Conductor limit
      end if

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

   end subroutine new_component_cpcm

end module moist_model_component_pcm_cpcm
