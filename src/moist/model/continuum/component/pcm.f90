!> Main PCM module - re-exports all PCM functionality
!>
!> - the top-level module to import for PCM access
module moist_model_continuum_component_pcm
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, &
      & pcm_solver_type, solver_type, moist_pcm_parameters_type
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, new_component_cpcm, &
      & model_continuum_component_cosmo, new_component_cosmo
   use moist_model_continuum_component_pcm_solvers, only: solve_pcm_lu, &
      & solve_pcm_cholesky, solve_pcm_iterative, solve_pcm_inversion, &
      & pcm_storage_auto, pcm_storage_dense, pcm_storage_matrix_free
   implicit none(type, external)
   public

end module moist_model_continuum_component_pcm
