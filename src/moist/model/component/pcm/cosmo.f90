!> COSMO (Conductor-like Screening Model) implementation
!> This module provides the COSMO variant of PCM with its specific dielectric
!> scaling (f epsilon = ( epsilon -1)/( epsilon +0.5))
module moist_model_component_pcm_cosmo
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_type, only: cavity_type, coupling_type
   use moist_model_component_pcm_type, only: pcm_base, solver_type, &
      & potential_source
   implicit none (type, external)
   private

   public :: cosmo
   public :: new_cosmo

   !> COSMO (Conductor-like Screening Model) variant
   !> Uses f epsilon = ( epsilon -1)/( epsilon +0.5) scaling
   type, extends(pcm_base) :: cosmo

   end type cosmo

contains

   !> Constructor for COSMO variant
   !> Sets f epsilon = ( epsilon -1)/( epsilon +0.5) and configures solver and potential source.
   subroutine new_cosmo(self, ctx, epsilon, solver, phi_source, &
                        & external_matrix, error)
      !> COSMO instance to initialize
      type(cosmo), intent(out) :: self
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
      ! negative, so the model is undefined there.
      if (epsilon < 1.0_wp) then
         call fatal_error(error, &
            & "[new_cosmo] Dielectric constant must be >= 1")
         return
      end if
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

      ! Handle external matrix
      if (present(external_matrix)) then
         call self%set_external_matrix(external_matrix)
      end if

      ! Set component name
      self%name = "COSMO"

   end subroutine new_cosmo

end module moist_model_component_pcm_cosmo
