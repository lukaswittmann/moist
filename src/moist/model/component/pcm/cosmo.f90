!> COSMO (Conductor-like Screening Model) implementation
!> This module provides the COSMO variant of PCM with its specific dielectric
!> scaling (f epsilon = ( epsilon -1)/( epsilon +0.5)). Matrix assembly is delegated to the cavity type.
module moist_model_component_pcm_cosmo
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_type, only: cavity_type, wavefunction_type
   use moist_model_component_pcm_type, only: pcm_base, solver_type, &
      & potential_source
   implicit none
   private

   public :: cosmo
   public :: new_cosmo

   !> COSMO (Conductor-like Screening Model) variant
   !> Uses f epsilon = ( epsilon -1)/( epsilon +0.5) scaling. Matrix assembly is delegated to the cavity type.
   type, extends(pcm_base) :: cosmo

      !> COSMO-specific: outlying charge correction factor (for D-COSMO)
      real(wp) :: outlying_charge_factor = 0.0_wp

   end type cosmo

contains

   !> Constructor for COSMO variant
   !> Sets f epsilon = ( epsilon -1)/( epsilon +0.5) and configures solver and potential source.
   subroutine new_cosmo(self, epsilon, solver, phi_source, &
                        & outlying_charge_factor, external_matrix, error)
      !> COSMO instance to initialize
      type(cosmo), intent(out) :: self
      !> Dielectric constant
      real(wp), intent(in) :: epsilon
      !> Optional: linear solver type
      integer, intent(in), optional :: solver
      !> Optional: potential source strategy
      integer, intent(in), optional :: phi_source
      !> Optional: outlying charge correction factor (D-COSMO)
      real(wp), intent(in), optional :: outlying_charge_factor
      !> Optional: external pre-computed matrix (ngrid, ngrid)
      real(wp), intent(in), optional :: external_matrix(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      ! Set dielectric properties
      self%epsilon = epsilon
      self%feps = (epsilon - 1.0_wp)/(epsilon + 0.5_wp)  ! COSMO formula

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

      ! Set COSMO-specific parameters
      if (present(outlying_charge_factor)) then
         self%outlying_charge_factor = outlying_charge_factor
      end if

      ! Handle external matrix
      if (present(external_matrix)) then
         call self%set_external_matrix(external_matrix)
      end if

      ! Set component name
      self%name = "COSMO"

   end subroutine new_cosmo

end module moist_model_component_pcm_cosmo
