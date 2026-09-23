!> COSMO (Conductor-like Screening Model) implementation
!>
!> The COSMO variant of PCM, with its specific dielectric scaling
!> (f epsilon = ( epsilon -1)/( epsilon +0.5))
module moist_model_component_pcm_cosmo
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_model_component_pcm_type, only: solvation_model_component_pcm, moist_pcm_parameters_type
   implicit none (type, external)
   private

   public :: solvation_model_component_cosmo
   public :: new_component_cosmo

   !> COSMO (Conductor-like Screening Model) variant
   !>
   !> - f epsilon = ( epsilon -1)/( epsilon +0.5) scaling
   type, extends(solvation_model_component_pcm) :: solvation_model_component_cosmo

   end type solvation_model_component_cosmo

contains

   !> Construct from parameter values; omission uses compiled defaults
   !>
   !> @param[inout] self Object to initialize
   !> @param[in] ctx Borrowed context; must outlive the object
   !> @param[in] epsilon Solvent dielectric constant
   !> @param[in] external_matrix Optional host-supplied PCM matrix
   !> @param[out] error Construction error
   !> @param[in] param Configuration copied by value
   subroutine new_component_cosmo(self, ctx, epsilon, external_matrix, error, param)
      !> COSMO instance to initialize
      type(solvation_model_component_cosmo), intent(out) :: self
      !> Shared run context (verbosity/debug/timer); borrowed, must outlive self
      type(moist_context_type), intent(in), target :: ctx
      !> Dielectric constant
      real(wp), intent(in) :: epsilon
      !> Solver configuration; omitted means compiled defaults
      type(moist_pcm_parameters_type), intent(in), optional :: param
      !> Optional: external pre-computed matrix (ngrid, ngrid)
      real(wp), intent(in), optional :: external_matrix(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Resolved solver configuration
      type(moist_pcm_parameters_type) :: settings

      if (present(param)) settings = param
      !> Borrow the shared run context (owns verbosity/debug/timer)
      self%ctx => ctx

      ! Set dielectric properties; below eps = 1 the scaling factor turns
      ! negative, so the model is undefined there
      if (epsilon < 1.0_wp) then
         call fatal_error(error, &
            & "[new_component_cosmo] Dielectric constant must be >= 1")
         return
      end if
      self%epsilon = epsilon
      self%feps = (epsilon - 1.0_wp)/(epsilon + 0.5_wp)  ! COSMO formula

      call settings%validate(error)
      if (allocated(error)) return
      self%solver = settings%solver
      self%solver_tol = settings%solver_tol
      self%solver_maxiter = settings%solver_maxiter

      ! Handle external matrix
      if (present(external_matrix)) then
         call self%set_external_matrix(external_matrix)
      end if

      ! Set component name
      self%name = "COSMO"

   end subroutine new_component_cosmo

end module moist_model_component_pcm_cosmo
